#!/usr/bin/env python3
"""Focused tests for the SdS areal-scan tail extractor."""

from pathlib import Path
import sys

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from analyze_sds_areal_scan import fit_tail_power, instantaneous_power, tail_model


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
    psi = np.column_stack([times**a, 2.0 * times**a, -times**a])
    pi = a * psi / times[:, None]
    diagnostic = instantaneous_power(times, psi, pi)
    np.testing.assert_allclose(diagnostic, a, rtol=2.0e-15, atol=2.0e-15)
    print("SdS areal-scan analysis tests passed")


if __name__ == "__main__":
    main()
