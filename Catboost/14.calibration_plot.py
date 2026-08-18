"""Draw one-vs-rest calibration curves from true labels and probabilities."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence
from SUtils import *
import numpy as np
from scipy import stats


def adjust_prob_logit(p, cutoff):
    """
    基于 best cutoff 的 logit 平移概率调整
    将 cutoff 映射到 0.5，保持 S 形曲线单调性
    """
    p = np.asarray(p, dtype=float)
    # 处理边界值避免除零
    p = np.clip(p, 1e-15, 1 - 1e-15)
    cutoff = np.clip(cutoff, 1e-15, 1 - 1e-15)

    # 闭式解，数值稳定
    p_new = (p * (1 - cutoff)) / (p * (1 - cutoff) + cutoff * (1 - p))
    return np.clip(p_new, 0, 1)


def _prepare_inputs(
    y_label: Sequence[object] | np.ndarray,
    y_pred: Sequence[float] | Sequence[Sequence[float]] | np.ndarray,
    class_names: Sequence[str] | None,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    labels = np.asarray(y_label)
    probabilities = np.asarray(y_pred, dtype=float)

    if labels.ndim != 1:
        raise ValueError("y_label must be a one-dimensional array.")
    if labels.size == 0:
        raise ValueError("y_label cannot be empty.")

    if probabilities.ndim == 2 and probabilities.shape[1] == 1:
        probabilities = probabilities[:, 0]
    if probabilities.ndim == 1:
        if probabilities.shape[0] != labels.shape[0]:
            raise ValueError("y_label and y_pred must contain the same number of samples.")
        probabilities = np.column_stack((1.0 - probabilities, probabilities))
    elif probabilities.ndim == 2:
        if probabilities.shape[0] != labels.shape[0]:
            raise ValueError("y_label and y_pred must contain the same number of samples.")
        if probabilities.shape[1] < 2:
            raise ValueError("y_pred must contain probabilities for at least two classes.")
    else:
        raise ValueError("y_pred must be a 1D positive-class probability or a 2D probability matrix.")

    if not np.isfinite(probabilities).all():
        raise ValueError("y_pred contains NaN or infinite values.")
    if ((probabilities < 0.0) | (probabilities > 1.0)).any():
        raise ValueError("All y_pred probabilities must be between 0 and 1.")
    if not np.allclose(probabilities.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError("Each row of a 2D y_pred matrix must sum to 1.")

    n_classes = probabilities.shape[1]
    if class_names is None:
        try:
            numeric_labels = labels.astype(float)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "class_names is required when y_label contains non-numeric labels."
            ) from exc
        if not np.equal(numeric_labels, np.floor(numeric_labels)).all():
            raise ValueError("Numeric y_label values must be integer class indices.")
        encoded_labels = numeric_labels.astype(int)
        classes = [str(index) for index in range(n_classes)]
    else:
        classes = [str(name) for name in class_names]
        if len(classes) != n_classes:
            raise ValueError(
                f"class_names has {len(classes)} entries, but y_pred has {n_classes} columns."
            )
        if len(set(classes)) != len(classes):
            raise ValueError("class_names must contain unique values.")

        label_strings = np.asarray([str(value) for value in labels], dtype=object)
        if set(label_strings).issubset(set(classes)):
            class_to_index = {name: index for index, name in enumerate(classes)}
            encoded_labels = np.asarray(
                [class_to_index[value] for value in label_strings], dtype=int
            )
        else:
            try:
                numeric_labels = labels.astype(float)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "Every y_label value must match class_names or be a zero-based class index."
                ) from exc
            if not np.equal(numeric_labels, np.floor(numeric_labels)).all():
                raise ValueError("Numeric y_label values must be integer class indices.")
            encoded_labels = numeric_labels.astype(int)

    if ((encoded_labels < 0) | (encoded_labels >= n_classes)).any():
        raise ValueError(f"y_label class indices must be between 0 and {n_classes - 1}.")

    return encoded_labels, probabilities, classes


def _calibration_table_from_prepared(
    encoded_labels: np.ndarray,
    probabilities: np.ndarray,
    classes: Sequence[str],
    n_bins: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    edges = np.linspace(0.0, 1.0, n_bins + 1)

    for class_index, class_name in enumerate(classes):
        class_probability = probabilities[:, class_index]
        bin_index = np.clip(
            np.digitize(class_probability, edges, right=True) - 1,
            0,
            n_bins - 1,
        )
        observed = (encoded_labels == class_index).astype(float)

        for current_bin in range(n_bins):
            mask = bin_index == current_bin
            if not mask.any():
                continue
            rows.append(
                {
                    "class_label": class_name,
                    "bin": current_bin + 1,
                    "n": int(mask.sum()),
                    "mean_predicted_probability": float(class_probability[mask].mean()),
                    "observed_fraction": float(observed[mask].mean()),
                }
            )

    return pd.DataFrame.from_records(rows)


def calibration_table(
    y_label: Sequence[object] | np.ndarray,
    y_pred: Sequence[float] | Sequence[Sequence[float]] | np.ndarray,
    class_names: Sequence[str] | None = None,
    n_bins: int = 10,
) -> pd.DataFrame:
    """Return the binned values used to draw one-vs-rest calibration curves.

    A one-dimensional ``y_pred`` is interpreted as the probability of
    ``class_names[1]`` in a binary problem. A two-dimensional ``y_pred`` must
    have one probability column per class, in the same order as ``class_names``.
    """
    if not isinstance(n_bins, int) or n_bins < 2:
        raise ValueError("n_bins must be an integer greater than or equal to 2.")
    encoded_labels, probabilities, classes = _prepare_inputs(
        y_label, y_pred, class_names
    )
    return _calibration_table_from_prepared(
        encoded_labels, probabilities, classes, n_bins
    )


def _validate_hl_inputs(
    y_true: Sequence[object] | np.ndarray,
    y_prob: Sequence[float] | np.ndarray,
    n_bins: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Validate and normalize the binary inputs used by the HL test."""
    if not isinstance(n_bins, int) or n_bins < 3:
        raise ValueError("n_bins must be an integer greater than or equal to 3.")

    labels = np.asarray(y_true)
    probabilities = np.asarray(y_prob, dtype=float)
    if labels.ndim != 1 or probabilities.ndim != 1:
        raise ValueError("y_true and y_prob must be one-dimensional arrays.")
    if labels.size == 0 or labels.size != probabilities.size:
        raise ValueError("y_true and y_prob must be non-empty and have the same length.")
    if not np.isfinite(probabilities).all() or ((probabilities < 0) | (probabilities > 1)).any():
        raise ValueError("y_prob must contain finite probabilities between 0 and 1.")

    try:
        label_values = labels.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("y_true must contain only binary values 0 and 1.") from exc
    if not np.isfinite(label_values).all() or not np.isin(label_values, (0, 1)).all():
        raise ValueError("y_true must contain only binary values 0 and 1.")
    labels = label_values.astype(int)
    if labels.size < n_bins:
        raise ValueError(f"Sample size ({labels.size}) must be >= n_bins ({n_bins}).")
    return labels, probabilities


def hosmer_lemeshow_test(
    y_true: Sequence[object] | np.ndarray,
    y_prob: Sequence[float] | np.ndarray,
    n_bins: int = 10,
) -> tuple[float, float]:
    """Return the standard Hosmer-Lemeshow chi-square statistic and P value.

    Observations are sorted by predicted probability and split into up to
    ``n_bins`` approximately equally sized groups. Each group contributes
    both the event and non-event terms to the statistic. The chi-square
    reference distribution uses ``valid_groups - 2`` degrees of freedom.
    """
    labels, probabilities = _validate_hl_inputs(y_true, y_prob, n_bins)
    order = np.argsort(probabilities, kind="mergesort")
    labels = labels[order]
    probabilities = probabilities[order]

    group_sizes = np.full(n_bins, labels.size // n_bins, dtype=int)
    group_sizes[: labels.size % n_bins] += 1

    chi2_stat = 0.0
    valid_groups = 0
    start = 0
    for group_size in group_sizes:
        stop = start + int(group_size)
        observed_events = float(labels[start:stop].sum())
        expected_events = float(probabilities[start:stop].sum())
        expected_non_events = float(group_size) - expected_events
        if expected_events <= 0.0 or expected_non_events <= 0.0:
            start = stop
            continue

        difference = observed_events - expected_events
        chi2_stat += difference * difference * float(group_size) / (
            expected_events * expected_non_events
        )
        valid_groups += 1
        start = stop

    if valid_groups < 3:
        return float("nan"), float("nan")

    p_value = float(stats.chi2.sf(chi2_stat, valid_groups - 2))
    return float(chi2_stat), p_value


def plot_calibration(
    y_label: Sequence[object] | np.ndarray,
    y_pred: Sequence[float] | Sequence[Sequence[float]] | np.ndarray,
    class_names: Sequence[str] | None = None,
    n_bins: int = 10,
    title: str = "OOF calibration",
    HL_chi2_stat: float | None = None,
    HL_p_value: float | None = None,
    output_png: str | Path | None = None,
    output_csv: str | Path | None = None,
    show: bool = False,
) -> tuple[plt.Figure, plt.Axes, pd.DataFrame]:
    """Create a calibration plot and run the Hosmer-Lemeshow test.

    For a binary model, the test uses the probability of class index 1. For a
    multiclass model, a separate one-vs-rest test is reported for every class.
    ``HL_chi2_stat`` and ``HL_p_value`` remain available for backward
    compatibility with callers that already calculated a binary result.
    """
    if not isinstance(n_bins, int) or n_bins < 3:
        raise ValueError("n_bins must be an integer greater than or equal to 3.")

    encoded_labels, probabilities, classes = _prepare_inputs(
        y_label, y_pred, class_names
    )
    table = _calibration_table_from_prepared(
        encoded_labels, probabilities, classes, n_bins
    )

    if (HL_chi2_stat is None) != (HL_p_value is None):
        raise ValueError("HL_chi2_stat and HL_p_value must be provided together.")

    hl_results: dict[str, tuple[float, float]] = {}
    if HL_chi2_stat is not None and HL_p_value is not None:
        if len(classes) != 2:
            raise ValueError("Precomputed HL values are only supported for binary models.")
        hl_results[classes[1]] = (float(HL_chi2_stat), float(HL_p_value))
    else:
        class_indices = [1] if len(classes) == 2 else range(len(classes))
        for class_index in class_indices:
            observed = (encoded_labels == class_index).astype(int)
            hl_results[classes[class_index]] = hosmer_lemeshow_test(
                observed, probabilities[:, class_index], n_bins
            )

    table.attrs["hosmer_lemeshow"] = {
        class_name: {"chi2_stat": chi2_stat, "p_value": p_value}
        for class_name, (chi2_stat, p_value) in hl_results.items()
    }
    for class_name, (chi2_stat, p_value) in hl_results.items():
        if len(hl_results) == 1:
            print(f"HL_chi2_stat: {chi2_stat}, HL_p_value: {p_value}")
        else:
            print(
                f"{class_name} - HL_chi2_stat: {chi2_stat}, "
                f"HL_p_value: {p_value}"
            )

    figure, axis = plt.subplots(figsize=(5, 4), constrained_layout=True)
    axis.plot([0, 1], [0, 1], "--", color="#666666", linewidth=1)
    for class_name, group in table.groupby("class_label", sort=False):
        axis.plot(
            group["mean_predicted_probability"],
            group["observed_fraction"],
            marker="o",
            label=str(class_name),
        )
    axis.set_xlabel("Mean predicted probability")
    axis.set_ylabel("Observed fraction")
    axis.set_title(title)
    axis.legend(loc="best")

    annotation_lines: list[str] = []
    for class_name, (chi2_stat, p_value) in hl_results.items():
        prefix = "" if len(hl_results) == 1 else f"{class_name}: "
        annotation_lines.append(f"{prefix}H-L chi2 = {chi2_stat:.3f}")
        annotation_lines.append(f"{prefix}H-L p = {p_value:.4f}")
    axis.text(
        0.97,
        0.03,
        "\n".join(annotation_lines),
        transform=axis.transAxes,
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="right",
    )

    if output_csv is not None:
        csv_path = Path(output_csv)
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(csv_path, index=False)
    if output_png is not None:
        png_path = Path(output_png)
        png_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(
            png_path,
            format="png",
            dpi=160,
            bbox_inches="tight",
            facecolor="white",
        )
    if show:
        plt.show()

    return figure, axis, table


if __name__ == '__main__':
    exp_names = ['subtype2_ALL_cluster1_vs_0', 'subtype2_ALL_cluster2_vs_0',
                'subtype2_Male_cluster1_vs_0', 'subtype2_Male_cluster2_vs_0',
                'subtype2_Female_cluster1_vs_0', 'subtype2_Female_cluster2_vs_0']
    model_names = ['Hae_Bio_PRS', 'Hae_Bio_PRS_covr10']
    n_bins=10
    for exp_Name in exp_names:
        for model_name in model_names:
            y_true, y_prob = [], []
            for fold in range(5):
                df_test_eval = pd.read_csv(multiview_result_path + "{}/{}/{}_test_eval.csv".format(fold, exp_Name, model_name))
                cutoff = get_best_cutoff(df_test_eval['label'].values, df_test_eval['pred_prob'].values)
                df_pos = df_test_eval[df_test_eval['label'] == 1]  # 正类（1）
                df_neg = df_test_eval[df_test_eval['label'] == 0]  # 负类（0）
                n_pos = len(df_pos)  # 正类样本数
                df_neg_sampled = df_neg.sample(n=n_pos, random_state=42)  # 从0类中随机抽取与1类等量的样本
                # 合并并重置索引（可选 shuffle 顺序）
                df_test_eval = pd.concat([df_pos, df_neg_sampled], ignore_index=True)
                print(df_test_eval.head(100))
                y_true_fold, y_prob_fold = df_test_eval['label'].values, df_test_eval['pred_prob'].values
                y_prob_fold = np.round(adjust_prob_logit(y_prob_fold, cutoff), 2)
                y_true.extend(y_true_fold)
                y_prob.extend(y_prob_fold)
            y_true, y_prob = np.array(y_true), np.array(y_prob)
            HL_chi2_stat, HL_p_value = hosmer_lemeshow_test(y_true, y_prob, n_bins)
            print(f"HL_chi2_stat: {HL_chi2_stat}, HL_p_value: {HL_p_value}")

            if y_prob.ndim == 1:
                # y_prob 形状为 (n_samples,)，假设是类别1（MDD）的概率
                y_prob = np.column_stack([1 - y_prob, y_prob])
            n_classes = y_prob.shape[1]

            plot_calibration(y_true, y_prob, ['HC', 'MDD'], n_bins,
            f"{exp_Name}_{model_name}", HL_chi2_stat, HL_p_value,
            os.path.join(multiview_result_path, exp_Name, f"{model_name}_Calibration.png"),
            os.path.join(multiview_result_path, exp_Name, f"{model_name}_Calibration.csv"),
            True)