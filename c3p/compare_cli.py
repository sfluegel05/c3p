from collections import defaultdict
from pathlib import Path
from typing import List, Optional, Annotated

import pandas as pd
import typer
from matplotlib import pyplot as plt
import seaborn as sns

from c3p.cli import configure_logging, verbose_option
from c3p.datamodel import Dataset, Config, ResultSet, EvaluationResult
from c3p.dumper import result_as_dict
from c3p.learn import safe_name
from c3p.outcome_stats import calculate_grouped_stats, assign_ranks_per_group

app = typer.Typer(help="Compare models.")

import logging

# Set up logger
logger = logging.getLogger()  # Get root logger

app = typer.Typer()


@app.command()
def combine(
        dirs: Optional[List[str]] = typer.Argument(..., help="directory name"),
        output_dir: Path = typer.Option(None, "--output", "-o", help="output directory"),
        exclude_incomplete: Optional[bool] = typer.Option(False, "--exclude-incomplete/--no-exclude-incomplete",
                                                          "-x/--no-x",
                                                          help="if set, exclude any result for a chemical not present in all sets"),
        verbose: Annotated[int, verbose_option] = 0
) -> None:
    """
    Evaluate a model on a dataset using a single class.
    """
    configure_logging(verbose)
    output_dir.mkdir(exist_ok=True, parents=True)
    rows = []
    for d in dirs:
        p = Path(d) / "eval"
        n = 0
        for fn in p.glob("*.json"):
            er = EvaluationResult.model_validate_json(fn.read_text())
            experiment_name = Path(d).name
            row = result_as_dict(er)
            row["experiment_name"] = experiment_name
            row["chemical_name"] = er.train_results.best_result.chemical_class.name
            rows.append(row)
            n += 1
        logger.info(f"Processed {n} results from {d}")
        if not n:
            logger.warning(f"No results found in {d}")
    if not rows:
        raise ValueError(f"No results found")
    df = pd.DataFrame(rows)
    expts = list(df["experiment_name"].unique())
    num_expts = len(expts)
    # if a chemical count is < num_expts, exclude it
    if exclude_incomplete:
        df = df.groupby("chemical_name").filter(lambda x: len(x) == num_expts)
    assign_ranks_per_group(df, "experiment_name")
    rank_stats = df.groupby('chemical_class').agg({
        'rank': ['mean', 'std', 'min', 'max']
    }).round(2)

    df_sorted = df.sort_values('rank')

    # Add range column
    rank_stats['rank', 'range'] = rank_stats['rank', 'max'] - rank_stats['rank', 'min']

    rank_stats_by_std = rank_stats.sort_values(('rank', 'std'))
    rank_stats_by_mean = rank_stats.sort_values(('rank', 'mean'))
    print("Rank stats by std")
    print(rank_stats_by_std)
    print("Rank stats by mean")
    print(rank_stats_by_mean)
    with open(output_dir / "rank_stats_by_std.csv", "w") as f:
        rank_stats_by_std.to_csv(f)
    with open(output_dir / "rank_stats_by_mean.csv", "w") as f:
        rank_stats_by_mean.to_csv(f)

    plt.figure(figsize=(8, 50))
    sns.boxplot(data=df_sorted, y='chemical_class', x='rank', order=rank_stats_by_mean.index)
    #plt.xticks(rotation=90, ha='right')
    plt.title('Rank Distribution by Chemical Class, sorted by mean rank')
    plt.tight_layout()
    plt.savefig(output_dir / "rank_distribution.png")

    top_n_variable = rank_stats.nsmallest(50, ('rank', 'mean')).index
    df_subset = df[df['chemical_class'].isin(top_n_variable)]
    plt.figure(figsize=(8, 12))
    sns.boxplot(data=df_subset, y='chemical_class', x='rank', order=top_n_variable)
    # plt.xticks(rotation=90, ha='right')
    plt.title('Rank Distribution, top 50')
    plt.tight_layout()
    plt.savefig(output_dir / "rank_distribution_top_50.png")

    plt.figure(figsize=(8, 50))
    sns.boxplot(data=df, y='chemical_class', x='rank', order=rank_stats_by_std.index)
    # plt.xticks(rotation=90, ha='right')
    plt.title('Rank Distribution by Chemical Class, sorted by std')
    plt.tight_layout()
    plt.savefig(output_dir / "rank_distribution_by_std.png")

    # Print summary stats
    print("Chemical classes with most consistent rankings (lowest std):")
    print(rank_stats.nsmallest(5, ('rank', 'std')))

    print("\nChemical classes with most variable rankings (highest std):")
    print(rank_stats.nlargest(5, ('rank', 'std')))

    # Optional: Calculate correlation between mean rank and rank variability
    correlation = rank_stats['rank']['mean'].corr(rank_stats['rank']['std'])
    print(f"\nCorrelation between mean rank and rank std: {correlation:.3f}")

    #rank_stats = rank_stats.reset_index()
    #rank_stats = rank_stats.sort_values(('std',))
    #print("Rank stats sorted by rank std")
    #print(rank_stats)
    print(df["experiment_name"].unique())
    stats = calculate_grouped_stats(df, "experiment_name")
    print(stats)
    print(stats.T)
    stats.to_csv(output_dir / "summary_stats.csv")
    stats.to_excel(output_dir / "summary_stats.xlsx")
    cols = ['micro_f1', 'macro_f1', 'micro_precision', 'macro_precision', 'micro_recall', 'macro_recall', 'micro_accuracy', 'macro_accuracy', 'micro_negative_predictive_value', 'macro_negative_predictive_value']
    stats_slim = stats[cols]
    stats_slim.to_csv(output_dir / "summary_stats_slim.csv")
    stats_slim.to_excel(output_dir / "summary_stats_slim.xlsx")

    stats.reset_index(inplace=True)
    stats.rename(columns={"index": "model"}, inplace=True)
    # rename values in model column, e.g. gpt-4o => 4o
    stats["model"] = stats["model"].str.replace("gpt-4o", "4o").str.replace("-preview", "").str.replace("-undef", "").str.replace("ensembl-5", "ensemble").str.replace("-hi", "-iter6").str.replace("-coder", "").str.replace("-3-sonnet", "")
    for typ in ["micro", "macro"]:
        main_metric = f"{typ}_f1"
        # Sort dataframe by 'f1' in descending order
        df_sorted = stats.sort_values(by=main_metric, ascending=False).reset_index(drop=True)

        # Plotting sorted bar chart
        fig, ax = plt.subplots(figsize=(10, 6))
        width = 0.25  # Bar width
        x = range(len(df_sorted))

        # Define distinct colors for the metrics
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']

        # Plot each metric as a separate group of bars with distinct colors
        ax.bar([i - width for i in x], df_sorted[main_metric], width, label=f'{typ} f1', color=colors[0])
        ax.bar(x, df_sorted[f'{typ}_precision'], width, label=f'{typ} precision', color=colors[1])
        ax.bar([i + width for i in x], df_sorted[f'{typ}_recall'], width, label=f'{typ} recall', color=colors[2])

        # Add labels and title
        ax.set_xlabel("Models")
        ax.set_ylabel("Scores")
        ax.set_title(f"Model Comparison for {typ} F1, Precision, and Recall (Sorted by {typ} F1)")
        ax.set_xticks(x)
        ax.set_xticklabels(df_sorted['model'])
        ax.legend()

        # Display the plot
        plt.tight_layout()
        #plt.show()
        plt.savefig(output_dir / f"model_comparison_{typ}.png")


if __name__ == "__main__":
    app()