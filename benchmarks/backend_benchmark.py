#!/usr/bin/env python3
"""Benchmark pyTALISES backend performance and parity metrics.

This script is designed to compare backend behavior on the same machine,
including GPU pods where both CPU and GPU backends are available.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import sys
import time
from dataclasses import asdict, dataclass

import numpy as np

# Support direct script execution from repo root:
#   python benchmarks/backend_benchmark.py ...
if __package__ in (None, ""):
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

import pytalises as pt


@dataclass
class BenchmarkResult:
    backend: str
    workload: str
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


def _make_problem(size: int):
    grid = pt.Grid(shape=(size,), extent=((-4.0, 4.0),))
    initial = ["exp(-x**2)", "0.2*exp(-(x-0.5)**2)"]
    potential = pt.HermitianPotential.from_lower_triangular(
        [
            "0.2*x**2 + 0.1*t",
            "0.03*cos(x) + 0.02*sin(t)",
            "0.1*x**2 - 0.05*t",
        ]
    )
    return grid, initial, potential


def _run_case(
    *,
    backend: str,
    workload: str,
    size: int,
    steps: int,
    dt: float,
    repeats: int,
    warmup: int,
) -> BenchmarkResult:
    timings: list[float] = []
    options = pt.PropagationOptions(backend=backend)

    for run_idx in range(warmup + repeats):
        grid, initial, potential = _make_problem(size)
        psi = pt.Wavefunction(initial, grid, normalize_const=1.0, backend=backend)

        _sync_backend(backend)
        t0 = time.perf_counter()
        if workload == "free":
            psi.freely_propagate(steps=steps, dt=dt, options=options)
        elif workload == "potential":
            psi.propagate(potential=potential, steps=steps, dt=dt, options=options)
        else:
            raise ValueError(f"Unknown workload '{workload}'")
        _sync_backend(backend)
        elapsed = time.perf_counter() - t0

        if run_idx < warmup:
            continue
        timings.append(elapsed)

    if not timings:
        raise RuntimeError("No benchmark timings collected; check repeats/warmup configuration.")

    arr = np.asarray(timings, dtype=float)
    return BenchmarkResult(
        backend=backend,
        workload=workload,
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

    grid, initial, potential = _make_problem(size)
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


def _collect_metadata() -> dict[str, object]:
    meta: dict[str, object] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "available_backends": list(pt.available_backends()),
    }

    if pt.has_cupy():
        import cupy as cp

        device = cp.cuda.Device()
        props = cp.cuda.runtime.getDeviceProperties(device.id)
        name_raw = props.get("name", b"unknown")
        if isinstance(name_raw, bytes):
            gpu_name = name_raw.decode("utf-8", errors="replace")
        else:
            gpu_name = str(name_raw)

        meta["cuda"] = {
            "device_id": int(device.id),
            "device_name": gpu_name,
            "cupy_version": cp.__version__,
        }

    return meta


def _parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark pyTALISES backends")
    parser.add_argument(
        "--sizes",
        default="192",
        help="Comma-separated 1D grid sizes (e.g. 128,192,256)",
    )
    parser.add_argument("--steps", type=int, default=25, help="Propagation steps")
    parser.add_argument("--dt", type=float, default=0.005, help="Time step")
    parser.add_argument("--repeats", type=int, default=3, help="Benchmark repeats")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup runs")
    parser.add_argument(
        "--workloads",
        default="free,potential",
        help="Comma-separated workloads: free,potential",
    )
    parser.add_argument("--json-out", type=Path, default=None, help="Optional JSON output path")
    parser.add_argument(
        "--min-speedup",
        type=float,
        default=None,
        help="Fail if CuPy speedup over NumPy is below this threshold",
    )
    args = parser.parse_args()

    if args.repeats < 1:
        parser.error("--repeats must be >= 1")
    if args.warmup < 0:
        parser.error("--warmup must be >= 0")

    sizes = [int(s) for s in _parse_csv(args.sizes)]
    workloads = _parse_csv(args.workloads)

    backends = ["numpy"]
    if pt.has_cupy():
        backends.append("cupy")

    results: list[BenchmarkResult] = []
    for workload in workloads:
        for size in sizes:
            for backend in backends:
                results.append(
                    _run_case(
                        backend=backend,
                        workload=workload,
                        size=size,
                        steps=args.steps,
                        dt=args.dt,
                        repeats=args.repeats,
                        warmup=args.warmup,
                    )
                )

    payload: dict[str, object] = {
        "metadata": _collect_metadata(),
        "results": [asdict(r) for r in results],
        "parity": _parity_check(
            size=max(sizes),
            steps=max(8, args.steps // 2),
            dt=args.dt,
        ),
    }

    print("Backend benchmark results:")
    for r in results:
        print(
            f"- {r.backend:>5} | {r.workload:>9} | n={r.size:>4}: "
            f"{r.mean_seconds:.4f}s ± {r.std_seconds:.4f}s"
        )

    speedups: dict[str, float] = {}
    if len(backends) >= 2:
        for workload in workloads:
            for size in sizes:
                numpy_row = next(
                    (
                        r
                        for r in results
                        if r.backend == "numpy" and r.workload == workload and r.size == size
                    ),
                    None,
                )
                cupy_row = next(
                    (
                        r
                        for r in results
                        if r.backend == "cupy" and r.workload == workload and r.size == size
                    ),
                    None,
                )
                if numpy_row and cupy_row:
                    key = f"{workload}:n{size}"
                    speedup = numpy_row.mean_seconds / cupy_row.mean_seconds
                    speedups[key] = speedup
                    print(f"- speedup {key} (numpy/cupy): {speedup:.2f}x")

    payload["speedups_numpy_over_cupy"] = speedups

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

    if args.min_speedup is not None:
        if "cupy" not in backends:
            print(
                "ERROR: --min-speedup requires CuPy benchmark results, "
                "but CuPy is unavailable on this host."
            )
            return 1
        if not speedups:
            print("ERROR: --min-speedup requested but no NumPy/CuPy comparisons were produced.")
            return 1
        if min(speedups.values()) < args.min_speedup:
            worst = min(speedups, key=speedups.get)
            print(
                f"ERROR: speedup {worst}={speedups[worst]:.2f}x is below required "
                f"minimum {args.min_speedup:.2f}x"
            )
            return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
