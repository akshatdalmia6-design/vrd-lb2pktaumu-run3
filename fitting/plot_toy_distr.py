import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, multivariate_normal
from itertools import product
import sys
from scipy.optimize import curve_fit
from scipy.stats import norm
from pathlib import Path

if len(sys.argv) < 2:
    print(f"usage: python {Path(sys.argv[0]).name} <toys.json>")
    sys.exit(1)

# File names
file_path = sys.argv[1]

with open(file_path, 'r') as f:
    toys_data = json.load(f)

# 2. Extract the "value" for Nsig, Nmid, and Ncom
# This assumes the JSON structure is a list of toy fits
params_to_plot = ["Nsig", "Ncom", "Expcom"]
data_dict = {p: [] for p in params_to_plot}

for toyid, toy in toys_data["toys"].items():
    for p in params_to_plot:
        if p in toy and "value" in toy[p]:
            data_dict[p].append(toy[p]["value"])

inputs = toys_data.get("inputs", {})


# 3. Define the Gaussian function for fitting
def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean)**2) / (2 * stddev**2))


# 4. Create the plots
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for i, p in enumerate(params_to_plot):
    data = np.array(data_dict[p], dtype=float)
    ax = axes[i]

    if data.size == 0:
        ax.set_title(f"{p}: no data")
        continue

    # Create histogram
    counts, bin_edges, _ = ax.hist(
        data,
        bins=30,
        color='skyblue',
        edgecolor='black',
        alpha=0.7,
        label='Toys')

    # Calculate bin centers for fitting
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    raw_mean, raw_std = float(np.mean(data)), float(np.std(data))
    textlines = [
        f"Sample mean: {raw_mean:.4g}", f"Sample std:  {raw_std:.4g}",
        f"N toys: {data.size}"
    ]

    # Initial guesses for the fit: [amplitude, mean, stddev]
    initial_guess = [max(counts), np.mean(data), np.std(data)]

    try:
        # Perform the fit
        popt, _ = curve_fit(
            gaussian, bin_centers, counts, p0=initial_guess, maxfev=10000)
        amp_fit, mean_fit, std_fit = popt

        # Plot the fitted curve
        x_fit = np.linspace(min(data), max(data), 100)
        y_fit = gaussian(x_fit, *popt)
        ax.plot(x_fit, y_fit, color='red', lw=2, label='Gaussian Fit')

        textlines += [
            f"Fit mean: {mean_fit:.4g}", f"Fit std:  {abs(std_fit):.4g}"
        ]

        print(
            f"Fit for {p}: Mean = {mean_fit:.4g}, StdDev = {abs(std_fit):.4g}")

    except Exception as e:
        print(f"Could not fit Gaussian for {p}: {e}")
        textlines.append("(Gaussian fit failed)")

    truth_key = {"Nsig": "Nsig", "Ncom": "Ncom", "Expcom": "Expcom"}[p]
    if truth_key in inputs and inputs[truth_key] is not None:
        ax.axvline(
            inputs[truth_key],
            color='black',
            ls='--',
            lw=2,
            label=f'injected = {inputs[truth_key]:.4g}')

    ax.text(
        0.05,
        0.95,
        "\n".join(textlines),
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.6))

    ax.set_title(f"Distribution of {p}")
    ax.set_xlabel("Value")
    ax.set_ylabel("Frequency")
    ax.legend()

output = Path(file_path).stem

plt.tight_layout()
plt.savefig(f'./toy_plots/{output}.pdf')
