#!/usr/bin/env python3
"""Plot and validate two ARCiS correlated-k tables at a chosen P and T."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits


def wavelength_edges(centres: np.ndarray) -> np.ndarray:
    centres = np.asarray(centres, dtype=np.float64)
    edges = np.empty(centres.size + 1, dtype=np.float64)
    edges[1:-1] = np.sqrt(centres[:-1] * centres[1:])
    edges[0] = centres[0] ** 2 / edges[1]
    edges[-1] = centres[-1] ** 2 / edges[-2]
    return edges


def gauss_legendre_weights(n: int) -> np.ndarray:
    _, weights = np.polynomial.legendre.leggauss(n)
    return np.asarray(weights / 2.0, dtype=np.float64)


def bracket_log(grid: np.ndarray, value: float) -> tuple[int, float, float]:
    """Return the bracketing index and ARCiS logarithmic grid weights."""
    grid = np.asarray(grid, dtype=np.float64)
    if grid.ndim != 1 or grid.size < 2 or np.any(np.diff(grid) <= 0.0):
        raise ValueError("Pressure and temperature grids must be increasing")
    if value <= 0.0:
        raise ValueError("Pressure and temperature must be positive")

    if value <= grid[0]:
        return 0, 1.0, 0.0
    if value >= grid[-1]:
        return grid.size - 2, 0.0, 1.0

    lower = int(np.searchsorted(grid, value, side="right")) - 1
    upper_weight = np.log10(value / grid[lower]) / np.log10(
        grid[lower + 1] / grid[lower]
    )
    return lower, 1.0 - upper_weight, upper_weight


def read_interpolated_table(
    path: Path, pressure: float, temperature: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """Read k(g, wavelength), using the interpolation in ReadOpacityFITS."""
    with fits.open(
        path,
        mode="readonly",
        memmap=True,
        do_not_scale_image_data=True,
    ) as hdul:
        data = hdul[0].data
        if data is None or data.ndim != 4:
            raise ValueError(f"{path}: HDU 0 is not a four-dimensional table")

        temperatures = np.asarray(hdul[1].data, dtype=np.float64).reshape(-1)
        pressures = np.asarray(hdul[2].data, dtype=np.float64).reshape(-1)
        wavelengths = np.asarray(hdul[3].data, dtype=np.float64).reshape(-1)

        i_pressure, wp0, wp1 = bracket_log(pressures, pressure)
        i_temperature, wt0, wt1 = bracket_log(temperatures, temperature)

        # FITS axes appear in reverse order in NumPy:
        # data[pressure, temperature, g, wavelength].
        ktable = (
            np.asarray(
                data[i_pressure, i_temperature, :, :], dtype=np.float64
            )
            * wp0
            * wt0
            + np.asarray(
                data[i_pressure, i_temperature + 1, :, :],
                dtype=np.float64,
            )
            * wp0
            * wt1
            + np.asarray(
                data[i_pressure + 1, i_temperature, :, :],
                dtype=np.float64,
            )
            * wp1
            * wt0
            + np.asarray(
                data[i_pressure + 1, i_temperature + 1, :, :],
                dtype=np.float64,
            )
            * wp1
            * wt1
        )

        if len(hdul) > 4 and hdul[4].data is not None:
            g_weights = np.asarray(hdul[4].data, dtype=np.float64).reshape(-1)
        else:
            g_weights = gauss_legendre_weights(ktable.shape[0])
        if g_weights.size != ktable.shape[0]:
            raise ValueError(f"{path}: the number of g weights is inconsistent")
        g_weights = g_weights / np.sum(g_weights)

        information = {
            "p0": pressures[i_pressure],
            "p1": pressures[i_pressure + 1],
            "wp0": wp0,
            "wp1": wp1,
            "t0": temperatures[i_temperature],
            "t1": temperatures[i_temperature + 1],
            "wt0": wt0,
            "wt1": wt1,
        }

    return wavelengths, ktable, g_weights, information


def wavelength_weighted_rebin(
    source_centres: np.ndarray,
    source_values: np.ndarray,
    target_centres: np.ndarray,
) -> np.ndarray:
    """Average scalar source values over the target wavelength cells."""
    source_edges = wavelength_edges(source_centres)
    target_edges = wavelength_edges(target_centres)
    result = np.full(target_centres.size, np.nan, dtype=np.float64)

    for index, (lower, upper) in enumerate(
        zip(target_edges[:-1], target_edges[1:])
    ):
        lower = max(lower, source_edges[0])
        upper = min(upper, source_edges[-1])
        if upper <= lower:
            continue

        first = max(
            0,
            int(np.searchsorted(source_edges, lower, side="right")) - 1,
        )
        stop = min(
            source_centres.size,
            int(np.searchsorted(source_edges, upper, side="left")),
        )
        cell_lower = source_edges[first:stop]
        cell_upper = source_edges[first + 1 : stop + 1]
        overlap = np.minimum(cell_upper, upper) - np.maximum(cell_lower, lower)
        keep = overlap > 0.0
        if np.any(keep):
            values = source_values[first:stop][keep]
            weights = overlap[keep]
            result[index] = np.sum(values * weights) / np.sum(weights)

    return result


def describe_interpolation(label: str, info: dict[str, float]) -> None:
    print(
        f"{label}: P bracket [{info['p0']:.6g}, {info['p1']:.6g}] bar "
        f"with weights [{info['wp0']:.5f}, {info['wp1']:.5f}]"
    )
    print(
        f"{label}: T bracket [{info['t0']:.6g}, {info['t1']:.6g}] K "
        f"with weights [{info['wt0']:.5f}, {info['wt1']:.5f}]"
    )


def compare(
    original_path: Path,
    reduced_path: Path,
    pressure: float,
    temperature: float,
    output_path: Path,
    show: bool,
) -> None:
    original_lam, original_k, original_wg, original_info = (
        read_interpolated_table(original_path, pressure, temperature)
    )
    reduced_lam, reduced_k, reduced_wg, reduced_info = (
        read_interpolated_table(reduced_path, pressure, temperature)
    )

    describe_interpolation("Original", original_info)
    describe_interpolation("Reduced ", reduced_info)

    original_mean = np.sum(original_k * original_wg[:, None], axis=0)
    reduced_mean = np.sum(reduced_k * reduced_wg[:, None], axis=0)
    original_mean_rebinned = wavelength_weighted_rebin(
        original_lam, original_mean, reduced_lam
    )

    valid = (
        np.isfinite(original_mean_rebinned)
        & np.isfinite(reduced_mean)
        & (original_mean_rebinned > 0.0)
    )
    relative_difference = np.full(reduced_mean.size, np.nan)
    relative_difference[valid] = (
        reduced_mean[valid] / original_mean_rebinned[valid] - 1.0
    )

    if np.any(valid):
        absolute_percent = 100.0 * np.abs(relative_difference[valid])
        print(
            "Weighted-mean comparison: "
            f"median |difference|={np.median(absolute_percent):.4g}%, "
            f"maximum |difference|={np.max(absolute_percent):.4g}%"
        )

    # ARCiS wavelength grids are stored in cm; convert to microns.
    original_micron = original_lam * 1.0e4
    reduced_micron = reduced_lam * 1.0e4

    import matplotlib.pyplot as plt

    figure, (opacity_axis, difference_axis) = plt.subplots(
        2,
        1,
        figsize=(11, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
        constrained_layout=True,
    )

    opacity_axis.loglog(
        original_micron,
        original_mean,
        color="0.65",
        linewidth=0.5,
        label="Original weighted mean",
    )
    opacity_axis.loglog(
        reduced_micron,
        original_mean_rebinned,
        color="black",
        linewidth=1.0,
        label="Original mean rebinned to reduced grid",
    )
    opacity_axis.loglog(
        reduced_micron,
        reduced_mean,
        color="tab:red",
        linewidth=0.9,
        linestyle="--",
        label="Saved reduced table",
    )
    opacity_axis.set_ylabel(r"Weighted mean $k$")
    opacity_axis.set_title(
        rf"ARCiS correlated-k comparison: $T={temperature:g}$ K, "
        rf"$P={pressure:g}$ bar"
    )
    opacity_axis.legend(loc="best")
    opacity_axis.grid(True, which="both", alpha=0.2)

    difference_axis.semilogx(
        reduced_micron,
        100.0 * relative_difference,
        color="tab:blue",
        linewidth=0.8,
    )
    difference_axis.axhline(0.0, color="black", linewidth=0.7)
    difference_axis.set_xlabel(r"Wavelength [$\mu$m]")
    difference_axis.set_ylabel("Difference [%]")
    difference_axis.grid(True, which="both", alpha=0.2)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180)
    print(f"Wrote {output_path}")
    if show:
        plt.show()
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare original and reduced ARCiS k-tables at one P,T"
    )
    parser.add_argument("original", type=Path)
    parser.add_argument("reduced", type=Path)
    parser.add_argument("--temperature", "-T", type=float, required=True)
    parser.add_argument("--pressure", "-P", type=float, required=True)
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("ktable_comparison.png"),
    )
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    compare(
        arguments.original,
        arguments.reduced,
        arguments.pressure,
        arguments.temperature,
        arguments.output,
        arguments.show,
    )
