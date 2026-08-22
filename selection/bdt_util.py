import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, roc_auc_score
from scipy import stats

DEFAULT_XGB_PARAMS = dict(
    n_estimators=2000,
    learning_rate=0.05,
    max_depth=4,
    min_child_weight=5.0,
    gamma=0.5,
    subsample=0.8,
    colsample_bytree=0.8,
    reg_lambda=5.0,
    reg_alpha=0.0,
    objective="binary:logistic",
    eval_metric="auc",
    early_stopping_rounds=80,
    importance_type="gain",
    tree_method="hist",
    n_jobs=4,
)


# ======================================================================
# Training
# ======================================================================
def train_kfold(X,
                y,
                weights,
                n_splits=5,
                params=None,
                strat_key=None,
                seed=42):

    params = {**DEFAULT_XGB_PARAMS, **(params or {})}
    X = np.asarray(X, dtype=np.float32)
    y = np.asarray(y)
    weights = np.asarray(weights, dtype=np.float32)
    strat = y if strat_key is None else np.asarray(strat_key)

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    oof_pred = np.full(len(y), np.nan)
    fold_id = np.full(len(y), -1, dtype=int)
    models = []

    for k, (tr, va) in enumerate(skf.split(X, strat)):
        model = xgb.XGBClassifier(**params, random_state=seed + k)
        model.fit(
            X[tr],
            y[tr],
            sample_weight=weights[tr],
            eval_set=[(X[va], y[va])],
            sample_weight_eval_set=[weights[va]],
            verbose=False,
        )
        best = model.best_iteration + 1
        oof_pred[va] = model.predict_proba(
            X[va], iteration_range=(0, best))[:, 1]
        fold_id[va] = k
        models.append(model)
        print(f"[fold {k}] best_iter={model.best_iteration:4d}  "
              f"val_AUC={model.best_score:.5f}")

    return models, oof_pred, fold_id


def save_models(models, path_stem):
    """Save each fold model as {path_stem}_fold{k}.json."""
    for k, m in enumerate(models):
        m.save_model(f"{path_stem}_fold{k}.json")
    print(f"saved {len(models)} models -> {path_stem}_fold*.json")


def load_models(path_stem, n_splits, params=None):
    models = []
    params = {**DEFAULT_XGB_PARAMS, **(params or {})}
    for k in range(n_splits):
        m = xgb.XGBClassifier(**params)
        m.load_model(f"{path_stem}_fold{k}.json")
        models.append(m)
    return models


# ======================================================================
# Diagnostics
# ======================================================================
def plot_roc(y, oof_pred, weights, out_path, fold_id=None):
    """Out-of-fold ROC + per-fold AUC spread (the honest performance number)."""
    y = np.asarray(y)
    weights = np.asarray(weights)
    auc = roc_auc_score(y, oof_pred, sample_weight=weights)
    fpr, tpr, _ = roc_curve(y, oof_pred, sample_weight=weights)

    title = f"out-of-fold AUC = {auc:.4f}"
    if fold_id is not None:
        fold_id = np.asarray(fold_id)
        aucs = []
        for k in np.unique(fold_id[fold_id >= 0]):
            m = fold_id == k
            aucs.append(
                roc_auc_score(y[m], oof_pred[m], sample_weight=weights[m]))
        title += f"   (per-fold {np.mean(aucs):.4f} +/- {np.std(aucs):.4f})"

    plt.figure(figsize=(5, 5))
    plt.plot(tpr, 1 - fpr, lw=2, color="k")
    plt.plot([0, 1], [1, 0], ls=":", c="grey", lw=1)
    plt.xlabel("Signal efficiency")
    plt.ylabel("Background rejection")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.title(title, fontsize=9)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.clf()
    print(f"saved {out_path}   AUC={auc:.4f}")
    return auc


def plot_importance(models, feat_names, out_path):
    """Mean gain importance across folds."""
    imp = np.mean([m.feature_importances_ for m in models], axis=0)
    order = np.argsort(imp)
    plt.figure(figsize=(6, 0.45 * len(feat_names) + 1))
    plt.barh(np.array(feat_names)[order], imp[order], color="tab:purple")
    plt.xlabel("mean gain importance")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.clf()
    ranked = dict(zip(np.array(feat_names)[order][::-1], imp[order][::-1]))
    print("feature importance (gain):")
    for n, v in ranked.items():
        print(f"  {n:<20} {v:.4f}")
    return ranked


def correlation_matrix(frame, names, out_path, title=""):
    """Linear correlation matrix (%) of `names` columns of a pandas frame."""
    M = np.corrcoef(
        np.vstack([np.asarray(frame[n], dtype=float) for n in names]))
    n = len(names)
    fig, ax = plt.subplots(figsize=(0.7 * n + 2, 0.7 * n + 2))
    im = ax.imshow(M, vmin=-1, vmax=1, cmap="coolwarm")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(names, fontsize=7)
    for i in range(n):
        for j in range(n):
            ax.text(
                j,
                i,
                f"{100*M[i, j]:.0f}",
                ha="center",
                va="center",
                fontsize=6,
                color="white" if abs(M[i, j]) > 0.5 else "black")
    fig.colorbar(im, fraction=0.046, pad=0.04)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"saved {out_path}")
    return M


def response_by_source(score, label, source_code, out_path, bins=30):
    """Validation: background BDT response on SS vs OS sideband separately.
    They should look similar; if not, the two proxies disagree and the
    'both combined' choice needs a rethink."""
    score = np.asarray(score)
    label = np.asarray(label)
    source_code = np.asarray(source_code)
    edges = np.linspace(0, 1, bins + 1)

    plt.figure(figsize=(6, 5))
    plt.hist(
        score[label == 1],
        bins=edges,
        density=True,
        histtype="step",
        lw=2,
        color="k",
        label="signal (sim)")
    for code, col, lab in [(0, "tab:red", "SS data"),
                           (1, "tab:orange", "OS sideband")]:
        sel = (label == 0) & (source_code == code)
        if sel.sum() == 0:
            continue
        plt.hist(
            score[sel],
            bins=edges,
            density=True,
            histtype="step",
            lw=2,
            color=col,
            label=lab)
    plt.yscale("log")
    plt.xlabel("BDT response")
    plt.ylabel("a.u.")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.clf()
    print(f"saved {out_path}")


def bdt_mass_correlation(scored_file,
                         out_path,
                         tree="tree",
                         mass_col="Lb_DTF_M_aux",
                         score_col="bdt",
                         label_col="label",
                         nbins_x=40,
                         nbins_y=25):

    import ROOT
    rdf = ROOT.RDataFrame(tree, str(scored_file)).Filter(f"{label_col}==0")
    arr = rdf.AsNumpy([mass_col, score_col])
    mass = np.asarray(arr[mass_col], dtype=float)
    score = np.asarray(arr[score_col], dtype=float)

    good = np.isfinite(mass) & np.isfinite(score)
    mass, score = mass[good], score[good]
    if len(mass) < 2:
        print("bdt_mass_correlation: <2 background events, skipping")
        return float("nan")

    rho = (float(np.corrcoef(mass, score)[0, 1])
           if mass.std() > 0 and score.std() > 0 else float("nan"))

    xlo, xhi = np.percentile(mass, [0.5, 99.5])
    H, xe, ye = np.histogram2d(
        mass, score, bins=[nbins_x, nbins_y], range=[[xlo, xhi], [0.0, 1.0]])
    colsum = H.sum(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        Hn = np.where(colsum > 0, H / colsum, 0.0)
    Hn_plot = np.ma.masked_where(Hn.T == 0, Hn.T)  # mask empty columns

    fig, ax = plt.subplots(figsize=(6, 5))
    pcm = ax.pcolormesh(xe, ye, Hn_plot, cmap="viridis")
    fig.colorbar(pcm, ax=ax, label="fraction / mass column")
    ax.set_xlabel("Lb_DTF_M [MeV]")
    ax.set_ylabel("BDT response")
    ax.set_title(f"background  (Pearson r = {rho:+.3f})", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"saved {out_path}   background BDT-vs-mass r = {rho:+.3f}")
    return rho


def sculpting_check(mass_bkg,
                    score_bkg,
                    out_path,
                    weights_bkg=None,
                    cuts=(0.0, 0.5, 0.8, 0.9),
                    mass_label="Lb_DTF_M [MeV]",
                    bins=50):
    """Does cutting on the BDT distort the background mass shape?
    Plots the (normalised) background mass for a few BDT thresholds and
    prints the linear correlation between score and mass. A flat family of
    shapes -> no sculpting; fanning shapes -> the BDT is sculpting the mass."""
    mass_bkg = np.asarray(mass_bkg, dtype=float)
    score_bkg = np.asarray(score_bkg, dtype=float)
    w = (np.ones_like(mass_bkg) if weights_bkg is None else np.asarray(
        weights_bkg, dtype=float))

    # weighted Pearson correlation
    def wcorr(a, b, w):
        wa = np.average(a, weights=w)
        wb = np.average(b, weights=w)
        cov = np.average((a - wa) * (b - wb), weights=w)
        va = np.average((a - wa)**2, weights=w)
        vb = np.average((b - wb)**2, weights=w)
        return cov / np.sqrt(va * vb) if va > 0 and vb > 0 else 0.0

    rho = wcorr(score_bkg, mass_bkg, w)
    print(f"corr(BDT score, mass) on background = {rho:+.3f}  "
          f"(near 0 is good)")

    lo, hi = np.percentile(mass_bkg, [0.5, 99.5])
    edges = np.linspace(lo, hi, bins + 1)
    plt.figure(figsize=(6, 5))
    for c in cuts:
        sel = score_bkg > c
        if sel.sum() < 10:
            continue
        plt.hist(
            mass_bkg[sel],
            bins=edges,
            weights=w[sel],
            density=True,
            histtype="step",
            lw=2,
            label=f"BDT > {c:.2f}")
    plt.xlabel(mass_label)
    plt.ylabel("normalised background")
    plt.title(f"sculpting check  (corr = {rho:+.3f})", fontsize=9)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.clf()
    print(f"saved {out_path}")
    return rho
