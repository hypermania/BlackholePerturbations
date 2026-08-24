#!/usr/bin/env python3
"""Focused tests for the SdS areal-scan tail extractor."""

from pathlib import Path
import sys
import tempfile

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from analyze_sds_areal_scan import (
    fit_tail_power,
    instantaneous_power,
    load_run,
    tail_model,
)


def write_run_layout(directory: Path, values: np.ndarray) -> None:
    (directory / "COMPLETE").write_text("complete\n", encoding="utf-8")
    (directory / "metadata.txt").write_text("run_name test\n", encoding="utf-8")
    np.array([0.0, 1.0], dtype=np.float64).tofile(directory / "t_list.dat")
    values.astype(np.float64).tofile(directory / "psi_list.dat")
    np.array([0.0], dtype=np.float64).tofile(directory / "x_grid.dat")
    np.array([0.0], dtype=np.float64).tofile(
        directory / "t_list_snapshots.dat"
    )


def main() -> None:
    times = np.linspace(400.0, 900.0, 20001)
    expected = np.array([-2.75, 175.0, 63.0])
    powers = tail_model(times, *expected)
    fit = fit_tail_power(times, powers, "synthetic", 400.0, 900.0)
    np.testing.assert_allclose(
        [fit.a, fit.b, fit.c], expected, rtol=2.0e-9, atol=2.0e-9
    )
    assert fit.identifiable

    # For psi=t^a, Pi=d_t psi, so the direct diagnostic must recover a.
    a = -3.25
    psi = times[:, None] ** a
    pi = a * psi / times[:, None]
    diagnostic = instantaneous_power(times, psi, pi)
    np.testing.assert_allclose(diagnostic, a, rtol=2.0e-15, atol=2.0e-15)

    with tempfile.TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory)
        single = root / "single"
        single.mkdir()
        write_run_layout(single, np.array([[1.0, 2.0], [3.0, 4.0]]))
        _, _, loaded_psi, loaded_pi = load_run(single)
        np.testing.assert_array_equal(loaded_psi[:, 0], [1.0, 3.0])
        np.testing.assert_array_equal(loaded_pi[:, 0], [2.0, 4.0])

        legacy = root / "legacy"
        legacy.mkdir()
        write_run_layout(
            legacy,
            np.array(
                [
                    [10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
                    [11.0, 21.0, 31.0, 41.0, 51.0, 61.0],
                ]
            ),
        )
        _, _, loaded_psi, loaded_pi = load_run(legacy)
        np.testing.assert_array_equal(loaded_psi[:, 0], [20.0, 21.0])
        np.testing.assert_array_equal(loaded_pi[:, 0], [50.0, 51.0])
    print("SdS areal-scan analysis tests passed")


if __name__ == "__main__":
    main()
