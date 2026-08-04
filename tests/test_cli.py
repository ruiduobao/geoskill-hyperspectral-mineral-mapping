"""CLI tests for hyperspectral-mineral-mapping."""
import subprocess
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod, SCRIPT


def run_cli(args, timeout=90):
    return subprocess.run([sys.executable, SCRIPT] + args,
                          capture_output=True, text=True, timeout=timeout)


class TestCLIBasics:
    def test_help(self):
        r = run_cli(["--help"])
        assert r.returncode == 0
        assert "--bbox" in r.stdout
        assert "--synthetic" in r.stdout
        assert "--library" in r.stdout

    def test_version(self):
        r = run_cli(["--version"])
        assert r.returncode == 0
        assert "1.0.0" in r.stdout

    def test_no_args_exit_2(self):
        assert run_cli([]).returncode == 2

    def test_bad_library_exit_2(self):
        r = run_cli(["--bbox", "116", "39", "117", "40", "--library", "bogus"])
        assert r.returncode == 2

    def test_input_not_found_exit_2(self):
        assert run_cli(["--input", "no_such.tif"]).returncode == 2

    def test_synthetic_without_bbox_exit_2(self):
        assert run_cli(["--synthetic"]).returncode == 2

    def test_synthetic_runs(self):
        out = "./_test_cli_out"
        r = run_cli(["--bbox", "116", "39", "117", "40", "--synthetic",
                     "--output-dir", out, "--quiet"])
        assert r.returncode == 0, r.stderr
        assert os.path.exists(os.path.join(out, "mineral_map.tif"))
        assert os.path.exists(os.path.join(out, "output-manifest.json"))

    def test_aster_library(self):
        out = "./_test_cli_aster"
        r = run_cli(["--bbox", "116", "39", "117", "40", "--synthetic",
                     "--library", "aster", "--output-dir", out, "--quiet"])
        assert r.returncode == 0, r.stderr

    def test_continuum_removal_flag(self):
        out = "./_test_cli_cr"
        r = run_cli(["--bbox", "116", "39", "117", "40", "--synthetic",
                     "--continuum-removal", "--output-dir", out, "--quiet"])
        assert r.returncode == 0, r.stderr
