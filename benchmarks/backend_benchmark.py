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


def _run_profiled_propagation(
    *,
    psi,
    workload: str,
    potential,
    steps: int,
    dt: float,
    options,
) -> dict[str, dict[str, float | int]]:
    if workload == "free":
        propagator = pt.Propagator(
            psi,
            potential=pt.DiagonalPotential(["0"] * psi.num_int_dim),
            variables={},
            options=options,
        )
        for _ in range(steps):
            propagator.kinetic_prop(dt)
        return propagator.stage_timings()

    if workload == "potential":
        propagator = pt.Propagator(
            psi,
            potential=potential,
            variables={},
            options=options,
        )
        propagator.kinetic_prop(dt / 2)
        propagator.potential_prop(dt)
        for _ in range(steps - 1):
            propagator.kinetic_prop(dt)
            propagator.potential_prop(dt)
        propagator.kinetic_prop(dt / 2)
        return propagator.stage_timings()

    raise ValueError(f"Unknown workload '{workload}'")


def _summarize_stage_runs(
    stage_runs: list[dict[str, dict[str, float | int]]],
    case_mean_seconds: float,
) -> dict[str, dict[str, float]]:
    if not stage_runs:
        return {}

    stage_names = sorted({name for run in stage_runs for name in run})
    summary: dict[str, dict[str, float]] = {}

    for name in stage_names:
        seconds = np.asarray(
            [float(run.get(name, {}).get("seconds", 0.0)) for run in stage_runs],
            dtype=float,
        )
        calls = np.asarray(
            [float(run.get(name, {}).get("calls", 0.0)) for run in stage_runs],
            dtype=float,
        )

        mean_seconds = float(np.mean(seconds))
        summary[name] = {
            "mean_seconds": mean_seconds,
            "std_seconds": float(np.std(seconds)),
            "mean_calls": float(np.mean(calls)),
            "share_of_case_runtime": (
                mean_seconds / case_mean_seconds if case_mean_seconds > 0 else 0.0
            ),
        }

    return summary


def _run_case(
    *,
    backend: str,
    workload: str,
    size: int,
    steps: int,
    dt: float,
    repeats: int,
    warmup: int,
    coupled_2x2_mode: str,
    potential_precompute_mode: str,
) -> tuple[BenchmarkResult, dict[str, dict[str, float]]]:
    timings: list[float] = []
    stage_runs: list[dict[str, dict[str, float | int]]] = []
    options = pt.PropagationOptions(
        backend=backend,
        profile_stages=True,
        coupled_2x2_mode=coupled_2x2_mode,
        potential_precompute_mode=potential_precompute_mode,
    )

    for run_idx in range(warmup + repeats):
        grid, initial, potential = _make_problem(size)
        psi = pt.Wavefunction(initial, grid, normalize_const=1.0, backend=backend)

        _sync_backend(backend)
        t0 = time.perf_counter()
        stage_profile = _run_profiled_propagation(
            psi=psi,
            workload=workload,
            potential=potential,
            steps=steps,
            dt=dt,
            options=options,
        )
        _sync_backend(backend)
        elapsed = time.perf_counter() - t0

        if run_idx < warmup:
            continue
        timings.append(elapsed)
        stage_runs.append(stage_profile)

    if not timings:
        raise RuntimeError("No benchmark timings collected; check repeats/warmup configuration.")

    arr = np.asarray(timings, dtype=float)
    result = BenchmarkResult(
        backend=backend,
        workload=workload,
        size=size,
        steps=steps,
        dt=dt,
        repeats=repeats,
        mean_seconds=float(np.mean(arr)),
        std_seconds=float(np.std(arr)),
    )
    stage_summary = _summarize_stage_runs(stage_runs, case_mean_seconds=result.mean_seconds)

    return result, stage_summary


def _parity_check(
    size: int,
    steps: int,
    dt: float,
    *,
    coupled_2x2_mode: str,
    potential_precompute_mode: str,
) -> dict[str, float] | None:
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
        options=pt.PropagationOptions(
            backend="numpy",
            coupled_2x2_mode=coupled_2x2_mode,
            potential_precompute_mode=potential_precompute_mode,
        ),
    )
    psi_cp.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(
            backend="cupy",
            coupled_2x2_mode=coupled_2x2_mode,
            potential_precompute_mode=potential_precompute_mode,
        ),
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
    parser.add_argument(
        "--coupled-2x2-mode",
        choices=("auto", "eigh"),
        default="auto",
        help="Potential-step strategy for 2x2 coupled systems",
    )
    parser.add_argument(
        "--potential-precompute-mode",
        choices=("off", "auto"),
        default="off",
        help="Time-affine potential precompute strategy",
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
    stage_breakdown: list[dict[str, object]] = []
    for workload in workloads:
        for size in sizes:
            for backend in backends:
                result, stages = _run_case(
                    backend=backend,
                    workload=workload,
                    size=size,
                    steps=args.steps,
                    dt=args.dt,
                    repeats=args.repeats,
                    warmup=args.warmup,
                    coupled_2x2_mode=args.coupled_2x2_mode,
                    potential_precompute_mode=args.potential_precompute_mode,
                )
                results.append(result)
                stage_breakdown.append(
                    {
                        "backend": backend,
                        "workload": workload,
                        "size": size,
                        "steps": args.steps,
                        "dt": args.dt,
                        "coupled_2x2_mode": args.coupled_2x2_mode,
                        "potential_precompute_mode": args.potential_precompute_mode,
                        "stages": stages,
                        "profiled_stage_seconds_total": float(
                            sum(stage["mean_seconds"] for stage in stages.values())
                        ),
                    }
                )

    payload: dict[str, object] = {
        "metadata": _collect_metadata(),
        "config": {
            "coupled_2x2_mode": args.coupled_2x2_mode,
            "potential_precompute_mode": args.potential_precompute_mode,
        },
        "results": [asdict(r) for r in results],
        "stage_breakdown": stage_breakdown,
        "parity": _parity_check(
            size=max(sizes),
            steps=max(8, args.steps // 2),
            dt=args.dt,
            coupled_2x2_mode=args.coupled_2x2_mode,
            potential_precompute_mode=args.potential_precompute_mode,
        ),
    }

    print("Backend benchmark results:")
    print(f"- coupled 2x2 mode: {args.coupled_2x2_mode}")
    print(f"- potential precompute mode: {args.potential_precompute_mode}")
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
