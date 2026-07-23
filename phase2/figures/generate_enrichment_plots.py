import os
import glob
import warnings

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

warnings.filterwarnings("ignore")

# -------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------
OUT_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/figures/"
DATA_PATTERN = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/results/aging_phase2.*.*.feature_enrichment.csv"

TOPIC_MAPS = {
    "Cell specificity": ["cell_specific"],
    "DNA Methylation Clock Sites": [
        "clocksites",
        "clocksites_Horvath2013",
        "clocksites_Shireby2020",
        "clocksites_Tong2024_BrainClock",
        "clocksites_Tong2024_Glia-In",
        "clocksites_Tong2024_Glia-Sin",
        "clocksites_Tong2024_Neu-In",
        "clocksites_Tong2024_Neu-Sin",
    ],
    "Feature size and complexity": ["gene_metrics", "peak_metrics"],
    "Genome Localization": ["telomere", "centromere"],
    "Human Accelerated Regions": ["har", "haqer"],
    "Regulatory Elements": [
        "Encode4_Proximal_enhancer",
        "Encode4_CA-CTCF",
        "Encode4_CA",
        "Encode4_Distal_enhancer",
        "Encode4_TF",
        "Encode4_CA-H3K4me34",
        "Encode4_Promoter",
        "TEs",
    ],
}

MIN_SIZE = 100
MAX_SIZE = 500


def get_source_to_topic_map():
    """Flattens the TOPIC_MAPS dictionary to map individual sources to their topics."""
    source_to_topic = {}
    for topic, sources in TOPIC_MAPS.items():
        for src in sources:
            source_to_topic[src] = topic
    return source_to_topic


def load_data(file_pattern, source_to_topic):
    """Loads all enrichment CSV files and combines them into a single DataFrame."""
    print("Loading data...")
    files = glob.glob(file_pattern)
    dfs = []

    for f in files:
        basename = os.path.basename(f)
        parts = basename.split(".")
        if len(parts) >= 4:
            modality = parts[1].upper()
            source = parts[2]

            if source in source_to_topic:
                try:
                    df = pd.read_csv(f)
                    df["modality"] = modality
                    df["source"] = source
                    df["Topic"] = source_to_topic[source]

                    if source in ["gene_metrics", "peak_metrics"]:
                        df = df[~df["test_name"].str.contains("chromosome::")]

                    dfs.append(df)
                except Exception as e:
                    print(f"Error reading {f}: {e}")

    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def clean_name(t):
    """Cleans up the test_name string for presentation in the plot."""
    t = str(t)
    t = t.replace("clock::", "").replace("encode::", "").replace("chromosome::", "Chr")
    t = t.replace("is_celltype_specific", "Cell-Type Specific")
    t = t.replace("clocksites_annotation_", "")
    t = t.replace("hg38_Encode4_cCRE.", "")
    t = t.replace("hg38_", "")
    return t


def process_data(df_all):
    """Filters data for significance, selects top features, and formats labels/metrics."""
    print("Filtering and processing metrics...")
    df_all["-log10p"] = -np.log10(df_all["p_value"].replace(0, 1e-300).astype(float))
    df_all["log2FC"] = np.log2(df_all["fold_change"].replace(0, 1e-10).astype(float))

    # Keep only significant hits
    df_sig = df_all[df_all["p_value"] < 0.05].copy()

    # Find the top 15 tests per topic by max -log10p
    test_scores = (
        df_sig.groupby(["Topic", "modality", "test_name"])["-log10p"]
        .max()
        .reset_index()
    )
    top_tests = []
    for topic in df_sig["Topic"].unique():
        top_t = test_scores[test_scores["Topic"] == topic].nlargest(15, "-log10p")
        top_tests.extend(top_t["test_name"].unique())

    df_plot = df_sig[df_sig["test_name"].isin(set(top_tests))].copy()

    df_plot["Clean_Name"] = df_plot["test_name"].apply(clean_name)
    df_plot["Y_label"] = (
        f"["
        + df_plot["Topic"]
        + "] "
        + df_plot["modality"]
        + " | "
        + df_plot["Clean_Name"]
    )

    # Bucket values for hue and size visual limits
    df_plot["log2FC_plot"] = df_plot["log2FC"].apply(
        lambda x: 3.5 if x > 3 else (-3.5 if x < -3 else x)
    )
    df_plot["log10p_plot"] = df_plot["-log10p"].apply(lambda x: 3.5 if x > 3 else x)

    # Sort to determine order
    df_plot.sort_values(
        by=["Topic", "modality", "Clean_Name"],
        ascending=[True, False, True],
        inplace=True,
    )

    return df_plot


def create_plot(df_plot, out_dir):
    """Generates and saves the refined enrichment bubble plot."""
    print("Generating Refined Bubble Plot...")

    y_order = df_plot["Y_label"].drop_duplicates().tolist()
    x_order = sorted(df_plot["cell_type"].unique())

    df_plot["cell_type"] = pd.Categorical(
        df_plot["cell_type"], categories=x_order, ordered=True
    )
    df_plot["Y_label"] = pd.Categorical(
        df_plot["Y_label"], categories=y_order, ordered=True
    )

    plt.figure(figsize=(18, max(12, len(y_order) * 0.4)))

    ax = sns.scatterplot(
        data=df_plot,
        x="cell_type",
        y="Y_label",
        size="log10p_plot",
        hue="log2FC_plot",
        palette="vlag",
        sizes=(MIN_SIZE, MAX_SIZE),
        edgecolor="gray",
        alpha=0.9,
        hue_norm=(-3.5, 3.5),
        size_norm=(0, 3.5),
        legend=False,
    )

    # Add horizontal separators between topics
    topics_in_order = [lbl.split("]")[0][1:] for lbl in y_order]
    last_topic = topics_in_order[0]
    for idx, current_topic in enumerate(topics_in_order):
        if current_topic != last_topic:
            plt.axhline(
                idx - 0.5, color="black", linewidth=1.5, linestyle="-", alpha=0.8
            )
            last_topic = current_topic

    plt.xticks(rotation=45, ha="right", fontsize=11)
    plt.yticks(fontsize=10)

    plt.title(
        "Multi-Modal Age-Associated Enrichments by Annotations (Significant Only)",
        fontsize=16,
        pad=20,
    )
    plt.xlabel("Cell Type", fontsize=14)
    plt.ylabel("Enrichment Feature", fontsize=14)

    plt.grid(True, linestyle=":", alpha=0.6, axis="x")
    plt.gca().invert_yaxis()

    # Build Custom Hue Legend
    cmap = sns.color_palette("vlag", as_cmap=True)
    norm_hue = plt.Normalize(-3.5, 3.5)

    hue_vals = [-3.5, -3, -2, -1, 0, 1, 2, 3, 3.5]
    hue_labels = ["< -3", "-3", "-2", "-1", "0", "1", "2", "3", "3+"]
    hue_handles = [
        mlines.Line2D(
            [],
            [],
            markerfacecolor=cmap(norm_hue(v)),
            markeredgecolor="gray",
            marker="o",
            linestyle="None",
            markersize=12,
        )
        for v in hue_vals
    ]

    # Build Custom Size Legend
    norm_size = plt.Normalize(0, 3.5)

    def get_markersize(v):
        scale = norm_size(v)
        area = MIN_SIZE + scale * (MAX_SIZE - MIN_SIZE)
        return np.sqrt(area)

    size_vals = [1, 2, 3, 3.5]
    size_labels = ["1", "2", "3", "3+"]
    size_handles = [
        mlines.Line2D(
            [],
            [],
            markerfacecolor="gray",
            markeredgecolor="gray",
            marker="o",
            linestyle="None",
            markersize=get_markersize(v),
            alpha=0.7,
        )
        for v in size_vals
    ]

    leg1 = ax.legend(
        hue_handles,
        hue_labels,
        title="Log2(Fold Change)",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        title_fontsize=12,
        fontsize=11,
        labelspacing=1.2,
    )
    ax.add_artist(leg1)

    ax.legend(
        size_handles,
        size_labels,
        title="-log10(p-value)",
        bbox_to_anchor=(1.02, 0.55),
        loc="upper left",
        frameon=False,
        title_fontsize=12,
        fontsize=11,
        labelspacing=1.5,
    )

    plt.tight_layout()

    out_png = os.path.join(out_dir, "enrichments_dot_plot.png")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")

    out_svg = os.path.join(out_dir, "enrichments_dot_plot.svg")
    plt.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close()

    print(f"Done! Saved outputs to:\n  - {out_png}\n  - {out_svg}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    source_to_topic = get_source_to_topic_map()
    df_all = load_data(DATA_PATTERN, source_to_topic)

    if df_all.empty:
        print("No data loaded. Please check the dataset path and pattern.")
        return

    df_plot = process_data(df_all)
    create_plot(df_plot, OUT_DIR)


if __name__ == "__main__":
    main()
