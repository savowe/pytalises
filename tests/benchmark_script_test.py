from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pytest

import benchmarks.backend_benchmark as bench


def test_benchmark_script_runs_directly_from_repo_root():
    repo_root = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [sys.executable, "benchmarks/backend_benchmark.py", "--help"],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0
    assert "Benchmark pyTALISES backends" in proc.stdout


def test_main_rejects_zero_repeats(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["backend_benchmark.py", "--repeats", "0"])
    with pytest.raises(SystemExit):
        bench.main()


def test_min_speedup_fails_when_cupy_unavailable(monkeypatch):
    monkeypatch.setattr(bench.pt, "has_cupy", lambda: False)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "backend_benchmark.py",
            "--sizes",
            "16",
            "--steps",
            "1",
            "--repeats",
            "1",
            "--warmup",
            "0",
            "--workloads",
            "free",
            "--min-speedup",
            "1.1",
        ],
    )

    assert bench.main() == 1
