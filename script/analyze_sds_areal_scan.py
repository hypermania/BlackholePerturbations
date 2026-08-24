#!/usr/bin/env python3
"""Analyze one sourced Schwarzschild-de Sitter areal-radius scan run."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar


OBSERVER_LABELS = ("x=50",)
TIME_SERIES_COLUMNS = ("psi_x50", "Pi_x50")
LEGACY_TIME_SERIES_COLUMN_COUNT = 6


@dataclass
class TailFit:
    observer: str
    fit_start: float
    fit_end: float
    a: float
    b: float
    c: float
    a_standard_error: float
    rms_residual: float
    point_count: int
    c_at_bound: bool
    jacobian_rank: int
    jacobian_condition_number: float
    identifiable: bool


@dataclass
class ExponentialDiagnostic:
    observer: str
    fit_start: float
    fit_end: float
    gamma: float
    r_squared: float
    point_count: int


def parse_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split(maxsplit=1)
        if fields:
            metadata[fields[0]] = fields[1] if len(fields) == 2 else ""
    return metadata


def load_run(run_directory: Path) -> tuple[
    dict[str, str], np.ndarray, np.ndarray, np.ndarray
]:
    required = (
        "COMPLETE",
        "metadata.txt",
        "t_list.dat",
        "psi_list.dat",
        "x_grid.dat",
        "t_list_snapshots.dat",
    )
    missing = [name for name in required if not (run_directory / name).is_file()]
    if missing:
        raise FileNotFoundError(
            f"incomplete run directory {run_directory}: missing {', '.join(missing)}"
        )

    metadata = parse_metadata(run_directory / "metadata.txt")
    times = np.fromfile(run_directory / "t_list.dat", dtype=np.float64)
    values = np.fromfile(run_directory / "psi_list.dat", dtype=np.float64)
    if times.size == 0:
        raise ValueError("time series is empty")
    if values.size == times.size * len(TIME_SERIES_COLUMNS):
        values = values.reshape(times.size, len(TIME_SERIES_COLUMNS))
        psi = values[:, :1]
        pi = values[:, 1:]
    elif values.size == times.size * LEGACY_TIME_SERIES_COLUMN_COUNT:
        # The completed pilot predates the single-observer output. Select its
        # x=50 columns so that the saved run remains reproducible.
        values = values.reshape(times.size, LEGACY_TIME_SERIES_COLUMN_COUNT)
        psi = values[:, 1:2]
        pi = values[:, 4:5]
    else:
        raise ValueError(
            "time-series size is inconsistent with the two-column output layout"
        )
    return metadata, times, psi, pi


def instantaneous_power(
    times: np.ndarray, psi: np.ndarray, pi: np.ndarray
) -> np.ndarray:
    result = np.full_like(psi, np.nan)
    valid = np.isfinite(psi) & np.isfinite(pi) & (np.abs(psi) > 1.0e-300)
    np.divide(times[:, None] * pi, psi, out=result, where=valid)
    return result


def tail_model(times: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    return a + b / (times - c)


def fit_tail_power(
    times: np.ndarray,
    powers: np.ndarray,
    observer: str,
    fit_start: float,
    fit_end: float,
) -> TailFit:
    mask = (
        (times >= fit_start)
        & (times <= fit_end)
        & np.isfinite(powers)
        & (np.abs(powers) < 1.0e4)
    )
    fit_times = times[mask]
    fit_powers = powers[mask]
    if fit_times.size < 100:
        raise ValueError(f"fewer than 100 valid fit points for {observer}")

    output_spacing = float(np.min(np.diff(fit_times)))
    # Dense output is useful for plots but unnecessary for a three-parameter fit.
    # Uniform thinning preserves the time weighting and reduces analysis time.
    if fit_times.size > 5000:
        selected = np.linspace(0, fit_times.size - 1, 5000).astype(int)
        fit_times = fit_times[selected]
        fit_powers = fit_powers[selected]

    span = fit_end - fit_start
    c_lower = -10.0 * max(fit_end, span, 1.0)
    c_upper = fit_start - max(0.5 * output_spacing, 1.0e-9)

    # For fixed c, a and b are linear parameters. Profile them analytically,
    # then minimize the one-dimensional residual in c. This remains reliable
    # when the full three-parameter Jacobian is nearly rank deficient.
    def profiled_fit(c_value: float) -> tuple[float, float, float]:
        inverse_time = 1.0 / (fit_times - c_value)
        centered_inverse = inverse_time - np.mean(inverse_time)
        denominator = float(centered_inverse @ centered_inverse)
        if denominator == 0:
            return float("inf"), float("nan"), float("nan")
        centered_power = fit_powers - np.mean(fit_powers)
        b_value = float(centered_inverse @ centered_power) / denominator
        a_value = float(np.mean(fit_powers) - b_value * np.mean(inverse_time))
        residual_value = a_value + b_value * inverse_time - fit_powers
        return float(residual_value @ residual_value), a_value, b_value

    minimum_distance = fit_start - c_upper
    maximum_distance = fit_start - c_lower
    distances = np.geomspace(minimum_distance, maximum_distance, 384)
    c_grid = fit_start - distances
    rss_grid = np.array([profiled_fit(c_value)[0] for c_value in c_grid])
    best_index = int(np.argmin(rss_grid))
    if best_index == 0 or best_index == c_grid.size - 1:
        c = float(c_grid[best_index])
    else:
        bracket_lower = float(min(c_grid[best_index - 1], c_grid[best_index + 1]))
        bracket_upper = float(max(c_grid[best_index - 1], c_grid[best_index + 1]))
        optimum = minimize_scalar(
            lambda c_value: profiled_fit(float(c_value))[0],
            bounds=(bracket_lower, bracket_upper),
            method="bounded",
            options={"xatol": 1.0e-12},
        )
        if not optimum.success:
            raise RuntimeError(f"tail fit failed for {observer}")
        c = float(optimum.x)
    _, a, b = profiled_fit(c)
    residual = tail_model(fit_times, a, b, c) - fit_powers
    jacobian = np.column_stack(
        (
            np.ones_like(fit_times),
            1.0 / (fit_times - c),
            b / (fit_times - c) ** 2,
        )
    )
    dof = max(fit_times.size - 3, 1)
    covariance = np.linalg.pinv(jacobian.T @ jacobian)
    covariance *= float(residual @ residual) / dof
    a_error = float(np.sqrt(max(covariance[0, 0], 0.0)))
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    jacobian_rank = int(np.linalg.matrix_rank(jacobian))
    condition_number = (
        float(singular_values[0] / singular_values[-1])
        if singular_values[-1] > 0
        else float("inf")
    )
    bound_tolerance = 1.0e-4 * max(abs(c_lower), abs(c_upper), 1.0)
    return TailFit(
        observer=observer,
        fit_start=fit_start,
        fit_end=fit_end,
        a=a,
        b=b,
        c=c,
        a_standard_error=a_error,
        rms_residual=float(np.sqrt(np.mean(residual**2))),
        point_count=int(fit_times.size),
        c_at_bound=(
            abs(c - c_lower) < bound_tolerance
            or abs(c - c_upper) < bound_tolerance
        ),
        jacobian_rank=jacobian_rank,
        jacobian_condition_number=condition_number,
        identifiable=(jacobian_rank == 3 and condition_number < 1.0e12),
    )


def run_fits(
    times: np.ndarray,
    powers: np.ndarray,
    fit_starts: tuple[float, ...],
    fit_end: float,
) -> list[TailFit]:
    return [
        fit_tail_power(times, powers[:, observer_index], observer, fit_start, fit_end)
        for fit_start in fit_starts
        for observer_index, observer in enumerate(OBSERVER_LABELS)
    ]


def fit_exponential_diagnostic(
    times: np.ndarray,
    powers: np.ndarray,
    fit_start: float = 250.0,
    fit_end: float = 450.0,
) -> list[ExponentialDiagnostic]:
    diagnostics = []
    for observer_index, observer in enumerate(OBSERVER_LABELS):
        absolute_power = np.abs(powers[:, observer_index])
        mask = (
            (times >= fit_start)
            & (times <= fit_end)
            & np.isfinite(absolute_power)
            & (absolute_power > 1.0e-300)
        )
        fit_times = times[mask]
        logarithm = np.log(absolute_power[mask])
        slope, intercept = np.polyfit(fit_times, logarithm, 1)
        predicted = intercept + slope * fit_times
        centered = logarithm - np.mean(logarithm)
        total_variation = float(centered @ centered)
        residual = logarithm - predicted
        r_squared = (
            1.0 - float(residual @ residual) / total_variation
            if total_variation > 0
            else float("nan")
        )
        diagnostics.append(
            ExponentialDiagnostic(
                observer=observer,
                fit_start=fit_start,
                fit_end=fit_end,
                gamma=float(-slope),
                r_squared=r_squared,
                point_count=int(fit_times.size),
            )
        )
    return diagnostics


def plot_loglog(
    output_path: Path,
    times: np.ndarray,
    psi: np.ndarray,
    metadata: dict[str, str],
    fit_start: float,
    fit_end: float,
) -> None:
    figure, axes = plt.subplots(
        1, len(OBSERVER_LABELS), figsize=(6, 4.5), sharex=True, squeeze=False
    )
    axes = axes[0]
    for observer_index, (axis, observer) in enumerate(zip(axes, OBSERVER_LABELS)):
        valid = (times > 0) & (np.abs(psi[:, observer_index]) > 0)
        axis.loglog(times[valid], np.abs(psi[valid, observer_index]), lw=0.8)
        axis.axvspan(fit_start, fit_end, color="tab:orange", alpha=0.12)
        axis.set_title(observer)
        axis.set_xlabel(r"$t$")
        axis.grid(alpha=0.25, which="both")
    axes[0].set_ylabel(r"$|\psi|$")
    figure.suptitle(
        rf"SdS areal source: $q={float(metadata['q_9LambdaM2']):g}$, "
        rf"$s={metadata['s']}$, $\ell={metadata['l']}$, "
        rf"$\beta={metadata['beta']}$"
    )
    figure.tight_layout()
    figure.savefig(output_path)
    plt.close(figure)


def plot_instantaneous_power(
    output_path: Path,
    times: np.ndarray,
    powers: np.ndarray,
    fits: list[TailFit],
    nominal_start: float,
    fit_end: float,
) -> None:
    figure, axes = plt.subplots(
        2, len(OBSERVER_LABELS), figsize=(7, 8), sharex="row", squeeze=False
    )
    for observer_index, observer in enumerate(OBSERVER_LABELS):
        overview_axis = axes[0, observer_index]
        fit_axis = axes[1, observer_index]
        valid = (
            (times >= 100)
            & np.isfinite(powers[:, observer_index])
            & (np.abs(powers[:, observer_index]) < 1.0e3)
        )
        overview_axis.plot(times[valid], powers[valid, observer_index], lw=0.7)
        fit = next(
            item
            for item in fits
            if item.observer == observer and item.fit_start == nominal_start
        )
        fit_times = np.linspace(nominal_start, fit_end, 800)
        fit_data = (
            (times >= nominal_start)
            & (times <= fit_end)
            & np.isfinite(powers[:, observer_index])
        )
        fit_axis.plot(
            times[fit_data], powers[fit_data, observer_index], lw=0.7,
            label=r"$p_{\rm loc}$",
        )
        fit_axis.plot(
            fit_times,
            tail_model(fit_times, fit.a, fit.b, fit.c),
            "--",
            lw=1.4,
            label=rf"fit $a={fit.a:.5g}$",
        )
        fit_axis.axhline(fit.a, color="black", ls=":", lw=1.0)
        fit_axis.set_yscale("symlog", linthresh=1.0e-27)
        overview_axis.axvspan(
            nominal_start, fit_end, color="tab:orange", alpha=0.12
        )
        overview_axis.set_title(observer)
        fit_axis.set_xlabel(r"$t$")
        overview_axis.grid(alpha=0.25)
        fit_axis.grid(alpha=0.25, which="both")
        fit_axis.legend()
        if not fit.identifiable:
            fit_axis.text(
                0.02,
                0.04,
                "three-parameter fit not identifiable",
                transform=fit_axis.transAxes,
                fontsize=8,
            )
    axes[0, 0].set_ylabel(r"$p_{\rm loc}=t\Pi/\psi$ (overview)")
    axes[1, 0].set_ylabel(r"$p_{\rm loc}$ (fit-window symlog)")
    figure.tight_layout()
    figure.savefig(output_path)
    plt.close(figure)


def plot_memory_residual(
    output_path: Path, times: np.ndarray, psi: np.ndarray
) -> np.ndarray:
    late_mask = times >= 0.9 * times[-1]
    memory = np.median(psi[late_mask], axis=0)
    residual = np.abs(psi - memory)
    figure, axes = plt.subplots(
        2, len(OBSERVER_LABELS), figsize=(7, 8), sharex=True, squeeze=False
    )
    for observer_index, observer in enumerate(OBSERVER_LABELS):
        valid = (times > 0) & (residual[:, observer_index] > 0)
        axes[0, observer_index].loglog(
            times[valid], residual[valid, observer_index], lw=0.7
        )
        axes[1, observer_index].semilogy(
            times[valid], residual[valid, observer_index], lw=0.7
        )
        axes[0, observer_index].set_title(
            rf"{observer}, $\psi_{{\rm mem}}={memory[observer_index]:.6g}$"
        )
        axes[0, observer_index].grid(alpha=0.25, which="both")
        axes[1, observer_index].grid(alpha=0.25, which="both")
        axes[1, observer_index].set_xlabel(r"$t$")
    axes[0, 0].set_ylabel(r"$|\psi-\psi_{\rm mem}|$ (log-log)")
    axes[1, 0].set_ylabel(r"$|\psi-\psi_{\rm mem}|$ (semilog)")
    figure.tight_layout()
    figure.savefig(output_path)
    plt.close(figure)
    return memory


def plot_snapshots(run_directory: Path, output_path: Path) -> None:
    x_grid = np.fromfile(run_directory / "x_grid.dat", dtype=np.float64)
    snapshot_times = np.fromfile(
        run_directory / "t_list_snapshots.dat", dtype=np.float64
    )
    if snapshot_times.size == 0:
        raise ValueError("no snapshots were saved")

    states = []
    for index in range(snapshot_times.size):
        state = np.fromfile(run_directory / f"state_{index}.dat", dtype=np.float64)
        if state.size != 2 * x_grid.size:
            raise ValueError(f"snapshot state_{index}.dat has the wrong size")
        states.append(state[: x_grid.size])
    states_array = np.asarray(states)

    figure, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    groups = (snapshot_times <= 200, snapshot_times >= 250)
    titles = ("Source passage and early response", "Intermediate and late response")
    for axis, group, title in zip(axes, groups, titles):
        indices = np.flatnonzero(group)
        colors = plt.cm.viridis(np.linspace(0.05, 0.95, max(indices.size, 1)))
        for color, index in zip(colors, indices):
            axis.plot(
                x_grid,
                states_array[index],
                color=color,
                lw=0.7,
                label=rf"$t={snapshot_times[index]:g}$",
            )
        axis.set_title(title)
        axis.set_ylabel(r"$\psi(t,x)$")
        axis.grid(alpha=0.2)
        axis.legend(ncol=5, fontsize=7, loc="best")
    axes[-1].set_xlabel(r"$x=r_*$")
    figure.tight_layout()
    figure.savefig(output_path)
    plt.close(figure)


def save_fit_results(
    run_directory: Path,
    metadata: dict[str, str],
    fits: list[TailFit],
    exponential_diagnostics: list[ExponentialDiagnostic],
    memory: np.ndarray,
    nominal_start: float,
) -> None:
    nominal = [fit for fit in fits if fit.fit_start == nominal_start]
    document = {
        "model": "p_loc(t) = a + b/(t-c)",
        "p_loc_definition": "d ln|psi| / d ln(t) = t Pi / psi",
        "tail_power_law": "a (signed; psi ~ t^a)",
        "run_name": metadata["run_name"],
        "memory_estimate_last_10_percent": {
            observer: float(value)
            for observer, value in zip(OBSERVER_LABELS, memory)
        },
        "nominal_fits": [asdict(fit) for fit in nominal],
        "window_sensitivity_fits": [asdict(fit) for fit in fits],
        "exponential_approach_diagnostic": [
            asdict(diagnostic) for diagnostic in exponential_diagnostics
        ],
    }
    (run_directory / "fit_results.json").write_text(
        json.dumps(document, indent=2) + "\n", encoding="utf-8"
    )

    lines = [
        "Model: p_loc(t) = a + b/(t-c)",
        "Definition: p_loc = d ln|psi| / d ln(t) = t Pi / psi",
        "The extracted signed tail power is a (psi ~ t^a).",
        "",
        "observer fit_start fit_end a a_std_error b c rms points "
        "c_at_bound jac_rank jac_condition identifiable",
    ]
    for fit in fits:
        lines.append(
            f"{fit.observer:>5s} {fit.fit_start:9.3f} {fit.fit_end:7.3f} "
            f"{fit.a: .12e} {fit.a_standard_error: .4e} "
            f"{fit.b: .12e} {fit.c: .12e} {fit.rms_residual: .4e} "
            f"{fit.point_count:d} {int(fit.c_at_bound):d} "
            f"{fit.jacobian_rank:d} {fit.jacobian_condition_number:.4e} "
            f"{int(fit.identifiable):d}"
        )
    lines.extend(
        [
            "",
            "Memory estimate (median over the final 10% of saved times):",
            *[
                f"{observer}: {value:.12e}"
                for observer, value in zip(OBSERVER_LABELS, memory)
            ],
            "",
            "Exponential approach diagnostic: |p_loc| ~ exp(-gamma t)",
            "observer fit_start fit_end gamma R_squared points",
            *[
                f"{item.observer:>5s} {item.fit_start:9.3f} "
                f"{item.fit_end:7.3f} {item.gamma:.12e} "
                f"{item.r_squared:.9f} {item.point_count:d}"
                for item in exponential_diagnostics
            ],
        ]
    )
    (run_directory / "fit_results.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def analyze_run(
    run_directory: Path,
    fit_starts: tuple[float, ...] = (400.0, 500.0, 600.0),
    fit_end: float = 900.0,
    nominal_start: float = 500.0,
) -> list[TailFit]:
    metadata, times, psi, pi = load_run(run_directory)
    if nominal_start not in fit_starts:
        raise ValueError("nominal fit start must be included in fit starts")
    if max(fit_starts) >= fit_end or fit_end > times[-1]:
        raise ValueError("fit windows do not lie inside the saved time range")

    powers = instantaneous_power(times, psi, pi)
    fits = run_fits(times, powers, fit_starts, fit_end)
    exponential_diagnostics = fit_exponential_diagnostic(times, powers)
    derived = np.column_stack((times, powers))
    derived.astype(np.float64).tofile(run_directory / "instantaneous_power.dat")
    (run_directory / "instantaneous_power_columns.txt").write_text(
        "t\np_loc_x50\n", encoding="utf-8"
    )

    plot_loglog(
        run_directory / "time_series_loglog.pdf",
        times,
        psi,
        metadata,
        nominal_start,
        fit_end,
    )
    plot_instantaneous_power(
        run_directory / "instantaneous_power.pdf",
        times,
        powers,
        fits,
        nominal_start,
        fit_end,
    )
    memory = plot_memory_residual(
        run_directory / "memory_subtracted.pdf", times, psi
    )
    plot_snapshots(run_directory, run_directory / "snapshots.pdf")
    save_fit_results(
        run_directory,
        metadata,
        fits,
        exponential_diagnostics,
        memory,
        nominal_start,
    )
    return fits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    parser.add_argument(
        "--fit-starts", type=float, nargs="+", default=(400.0, 500.0, 600.0)
    )
    parser.add_argument("--fit-end", type=float, default=900.0)
    parser.add_argument("--nominal-start", type=float, default=500.0)
    arguments = parser.parse_args()

    fits = analyze_run(
        arguments.run_directory,
        tuple(arguments.fit_starts),
        arguments.fit_end,
        arguments.nominal_start,
    )
    print(f"Analyzed {arguments.run_directory}")
    for fit in fits:
        if fit.fit_start == arguments.nominal_start:
            print(
                f"  {fit.observer}: a={fit.a:.12g}, "
                f"b={fit.b:.12g}, c={fit.c:.12g}, "
                f"rms={fit.rms_residual:.3g}"
            )


if __name__ == "__main__":
    main()
