import argparse
import glob
import json
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


def injected_from_name(path):
    """Pull the injected nsig out of the filename, for diagnostics only."""
    m = re.search(r"nsig(\d+(?:\.\d+)?)", path)
    return float(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pattern", help="glob for the FC JSONs (QUOTE it)")
    ap.add_argument(
        "--observed",
        type=float,
        default=None,
        help="observed N_sig_hat, drawn as a vertical line")
    ap.add_argument("--out", default="plots/empirical_belt.pdf")
    ap.add_argument(
        "--bin",
        type=float,
        default=0.0,
        help="if > 0, average runs into bins of this width in "
        "measured N_sig (smooths the edges)")
    args = ap.parse_args()

    files = sorted(glob.glob(args.pattern))
    if not files:
        raise SystemExit(f"no files matched {args.pattern!r}")

    rows, skipped = [], 0
    for f in files:
        d = json.load(open(f))
        if d.get("upp_nsig") is None:
            skipped += 1
            continue
        rows.append(
            (float(d["val_nsig"]),
             float(d["low_nsig"] if d["low_nsig"] is not None else 0.0),
             float(d["upp_nsig"]), injected_from_name(f)))

    if not rows:
        raise SystemExit("no run produced a crossing -- widen --fcbins")

    rows.sort(key=lambda r: r[0])  # sort by MEASURED value
    x = np.array([r[0] for r in rows])
    lo = np.array([r[1] for r in rows])
    hi = np.array([r[2] for r in rows])
    inj = [r[3] for r in rows]

    print(f"runs used {len(rows)}   no crossing {skipped}")
    print(f"\n{'injected':>9} {'measured':>9} {'low':>8} {'upper':>8}")
    for (xi, li, hi_, ij) in rows:
        tag = f"{ij:9.1f}" if ij is not None else "        ?"
        print(f"{tag} {xi:9.2f} {li:8.3f} {hi_:8.3f}")

    # #fix5-style smoothing: several runs can land on the same measured value
    # (the boundary piles them up at 0), and the edges are noisy because each
    # run is a different dataset with different nuisances.
    if args.bin > 0:
        edges = np.arange(x.min(), x.max() + args.bin, args.bin)
        idx = np.digitize(x, edges) - 1
        xb, lob, hib = [], [], []
        for k in np.unique(idx):
            m = idx == k
            xb.append(x[m].mean())
            lob.append(lo[m].mean())
            hib.append(hi[m].mean())
        x, lo, hi = np.array(xb), np.array(lob), np.array(hib)
        print(f"\nbinned into {len(x)} points of width {args.bin}")

    # monotonicity check: the belt edges should rise with the measurement
    if np.any(np.diff(hi) < 0):
        n_bad = int((np.diff(hi) < 0).sum())
        print(
            f"\nNOTE: the upper edge decreases at {n_bad} step(s). Each point "
            "is a different dataset, so some non-monotonicity is expected; "
            "many such steps means more runs (or --bin) are needed.")

    plt.figure(figsize=(6.5, 5.5))
    plt.fill_between(
        x, lo, hi, color="tab:blue", alpha=0.30, label="95% CL belt")
    plt.plot(x, hi, color="tab:blue", lw=2)
    plt.plot(x, lo, color="tab:blue", lw=2)

    span = [min(x.min(), lo.min()), max(x.max(), hi.max())]
    plt.plot(span, span, ls=":", color="grey", lw=1.5, label="measured = true")

    if args.observed is not None:
        plt.axvline(
            args.observed,
            color="red",
            lw=2,
            label=rf"observed $\hat N_{{\rm sig}}$ = {args.observed:.2f}")
        u = float(np.interp(args.observed, x, hi))
        plt.axhline(
            u, color="red", ls="--", lw=2, label=f"upper limit = {u:.2f}")
        print(f"\nread off at observed = {args.observed:.2f}: "
              f"upper limit = {u:.3f}")

    ax = plt.gca()
    ax.set_xlim(0.0, 20.0)
    ax.set_ylim(0.0, 20.0)
    ax.set_aspect("equal", adjustable="box")

    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(MultipleLocator(5.0))
        axis.set_minor_locator(MultipleLocator(1.0))
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.tick_params(which="major", length=7)
    ax.tick_params(which="minor", length=3.5)

    plt.xlabel(r"$\hat N_{\rm sig}$  (measured)")
    plt.ylabel(r"$N_{\rm sig}$  (true)")
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(args.out)
    print(f"\nCreated: {args.out}")


if __name__ == "__main__":
    main()
