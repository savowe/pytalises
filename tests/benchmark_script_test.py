from __future__ import annotations

import json
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


def test_stage_breakdown_is_emitted_in_json(monkeypatch, tmp_path):
    json_out = tmp_path / "benchmark.json"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "backend_benchmark.py",
            "--sizes",
            "16",
            "--steps",
            "2",
            "--repeats",
            "1",
            "--warmup",
            "0",
            "--workloads",
            "free",
            "--json-out",
            str(json_out),
        ],
    )

    assert bench.main() == 0

    payload = json.loads(json_out.read_text(encoding="utf-8"))
    breakdown = payload.get("stage_breakdown")
    assert isinstance(breakdown, list)
    assert breakdown

    entry = breakdown[0]
    assert entry["backend"] == "numpy"
    assert entry["workload"] == "free"
    assert "kinetic_prop" in entry["stages"]
    assert entry["stages"]["kinetic_prop"]["mean_seconds"] >= 0.0
