#!/usr/bin/env python3
"""Benchmark pyTALISES backend performance and basic parity metrics."""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

import pytalises as pt


@dataclass
class BenchmarkResult:
    backend: str
    size: int
    steps: int
    dt: float
    repeats: int
    mean_seconds: float
    std_seconds: float


def _sync_backend(backend: str) -> None:
    if backend != "cupy":
        return
    import cupy as cp

    cp.cuda.Stream.null.synchronize()


def _run_case(*, backend: str, size: int, steps: int, dt: float, repeats: int, warmup: int) -> BenchmarkResult:
    grid = pt.Grid(shape=(size,), extent=((-4.0, 4.0),))
    initial = ["exp(-x**2)", "0.2*exp(-(x-0.5)**2)"]
    potential = pt.HermitianPotential.from_lower_triangular(
        [
            "0.2*x**2 + 0.1*t",
            "0.03*cos(x) + 0.02*sin(t)",
            "0.1*x**2 - 0.05*t",
        ]
    )

    timings: list[float] = []
    options = pt.PropagationOptions(backend=backend)

    for run_idx in range(warmup + repeats):
        psi = pt.Wavefunction(initial, grid, normalize_const=1.0, backend=backend)
        _sync_backend(backend)
        t0 = time.perf_counter()
        psi.propagate(potential=potential, steps=steps, dt=dt, options=options)
        _sync_backend(backend)
        elapsed = time.perf_counter() - t0

        if run_idx < warmup:
            continue
        timings.append(elapsed)

    arr = np.asarray(timings, dtype=float)
    return BenchmarkResult(
        backend=backend,
        size=size,
        steps=steps,
        dt=dt,
        repeats=repeats,
        mean_seconds=float(np.mean(arr)),
        std_seconds=float(np.std(arr)),
    )


def _parity_check(size: int, steps: int, dt: float) -> dict[str, float] | None:
    if not pt.has_cupy():
        return None

    import cupy as cp

    grid = pt.Grid(shape=(size,), extent=((-4.0, 4.0),))
    initial = ["exp(-x**2)", "0.2*exp(-(x-0.5)**2)"]
    potential = pt.HermitianPotential.from_lower_triangular(
        [
            "0.2*x**2 + 0.1*t",
            "0.03*cos(x) + 0.02*sin(t)",
            "0.1*x**2 - 0.05*t",
        ]
    )

    psi_np = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="numpy")
    psi_cp = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="cupy")

    psi_np.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(backend="numpy"),
    )
    psi_cp.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(backend="cupy"),
    )

    amp_np = np.asarray(psi_np.amp)
    amp_cp = cp.asnumpy(psi_cp.amp)

    density_diff = np.abs(np.abs(amp_np) ** 2 - np.abs(amp_cp) ** 2)
    occ_np = np.asarray(psi_np.state_occupation())
    occ_cp = cp.asnumpy(psi_cp.state_occupation())

    return {
        "max_density_abs_diff": float(np.max(density_diff)),
        "max_occupation_abs_diff": float(np.max(np.abs(occ_np - occ_cp))),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark pyTALISES backends")
    parser.add_argument("--size", type=int, default=192, help="1D grid size")
    parser.add_argument("--steps", type=int, default=25, help="Propagation steps")
    parser.add_argument("--dt", type=float, default=0.005, help="Time step")
    parser.add_argument("--repeats", type=int, default=3, help="Benchmark repeats")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup runs")
    parser.add_argument("--json-out", type=Path, default=None, help="Optional JSON output path")
    parser.add_argument(
        "--min-speedup",
        type=float,
        default=None,
        help="Fail if CuPy speedup over NumPy is below this threshold",
    )
    args = parser.parse_args()

    results = [
        _run_case(
            backend="numpy",
            size=args.size,
            steps=args.steps,
            dt=args.dt,
            repeats=args.repeats,
            warmup=args.warmup,
        )
    ]

    if pt.has_cupy():
        results.append(
            _run_case(
                backend="cupy",
                size=args.size,
                steps=args.steps,
                dt=args.dt,
                repeats=args.repeats,
                warmup=args.warmup,
            )
        )

    payload: dict[str, object] = {
        "results": [asdict(r) for r in results],
        "parity": _parity_check(size=args.size, steps=max(8, args.steps // 2), dt=args.dt),
    }

    print("Backend benchmark results:")
    for r in results:
        print(f"- {r.backend:>5}: {r.mean_seconds:.4f}s ± {r.std_seconds:.4f}s")

    speedup = None
    if len(results) == 2:
        numpy_t = next(r.mean_seconds for r in results if r.backend == "numpy")
        cupy_t = next(r.mean_seconds for r in results if r.backend == "cupy")
        speedup = numpy_t / cupy_t
        payload["speedup_numpy_over_cupy"] = speedup
        print(f"- speedup (numpy/cupy): {speedup:.2f}x")

    parity = payload.get("parity")
    if isinstance(parity, dict):
        print(
            "- parity max abs diff: density={:.3e}, occupation={:.3e}".format(
                parity["max_density_abs_diff"],
                parity["max_occupation_abs_diff"],
            )
        )

    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    if args.min_speedup is not None and speedup is not None and speedup < args.min_speedup:
        print(
            f"ERROR: speedup {speedup:.2f}x is below required minimum {args.min_speedup:.2f}x"
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
