import math
import os
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Annotated

import pandas as pd
import typer

from c3p.cli import verbose_option, configure_logging
from c3p.datamodel import Dataset, Config, EvaluationResult, ResultSet, CodeStatistics
from c3p.dumper import result_as_dict
from c3p.learn import safe_name
from c3p.outcome_stats import calculate_all_outcome_stats
from c3p.plotting import plot_scatter

app = typer.Typer(help="CHEBI learn.")

import logging

# Set up logger
logger = logging.getLogger()  # Get root logger


app = typer.Typer()


@app.command()
def summarize(
        working_dir: Path = typer.Option(None, "--workdir", "-w", help="path to workdir"),
        output_dir: Optional[Path] = typer.Option(None, "--output", "-o", help="output directory"),
        verbose: Annotated[int, verbose_option] = 0
) -> None:
    """
    Evaluate a model on a dataset using a single class.
    """
    n = 0
    configure_logging(verbose)
    eval_results = []
    eval_objects = []
    for fn in working_dir.glob("*.json"):
        if "/cache" in str(working_dir):
            er = ResultSet.model_validate_json(fn.read_text())
        else:
            er = EvaluationResult.model_validate_json(fn.read_text())
        eval_objects.append(er)
        #auprc = calculate_auprc(er)
        er_dict = result_as_dict(er)
        if isinstance(er, EvaluationResult):
            br = er.train_results.best_result
            if not br and "chebifier" in str(working_dir):
                br = er.test_result
        else:
            br = er
        if not br.code_statistics and br.code:
            print(f"calculated code statistics for {fn}")
            br.code_statistics = CodeStatistics.from_code(br.code)
            for k, v in br.code_statistics.model_dump().items():
                if isinstance(v, (int, float)):
                    er_dict[k] = v
        er_dict["num_train_examples"] = br.num_true_positives + br.num_false_negatives
        #if not er_dict["num_train_examples"]:
        #    print(f"no train examples for {er_dict['f1']} {fn}")
        timestamp = os.path.getmtime(fn)
        # Convert to datetime object
        dt = datetime.fromtimestamp(timestamp)
        formatted_time = dt.strftime("%Y-%m-%d %H:%M:%S")
        er_dict["time"] = formatted_time
        #er_dict["auprc"] = auprc
        eval_results.append(er_dict)
    df = pd.DataFrame(eval_results)
    df["log_num_train_examples"] = df["num_train_examples"].apply(lambda x: 0 if x == 0 else math.log(x) / math.log(10))
    print(df.describe())
    metrics = ["f1", "precision", "recall", "accuracy", "negative_predictive_value"]
    print(df.aggregate({k: ["mean", "std", "max", "min"] for k in metrics}))
    all_stats_sf = calculate_all_outcome_stats(df)
    print(all_stats_sf.T.round(3))
    if output_dir:
        output_dir.mkdir(exist_ok=True, parents=True)
        md_dir = output_dir / "markdown"
        md_dir.mkdir(exist_ok=True, parents=True)
        for eo in eval_objects:
            br = eo.train_results.best_result
            md_fn = md_dir / f"{safe_name(br.chemical_class.name)}.md"
            with open(md_fn, "w") as f:
                f.write(eo.markdown)
        df.to_csv(output_dir / "results.csv", index=False)
        df.to_excel(output_dir / "results.xslx", index=False)
        all_stats_sf.to_csv(output_dir / "outcome_stats.csv", index=False)
        df.describe().to_csv(output_dir / "summary.csv", index=True)
        df.describe().to_html(output_dir / "summary.html", index=True)
        plot_scatter(df, "train_f1", "f1", results_dir=output_dir)
        plot_scatter(df, "log_num_train_examples", "f1", results_dir=output_dir)
        plot_scatter(df, "complexity", "f1", results_dir=output_dir)
        #plot_scatter(df, "num_train_examples", "f1", results_dir=output_dir)
        for min_f1 in [0.05, 0.5]:
            filtered_df = df[df["train_f1"] > min_f1]
            filtered_df.to_csv(output_dir / f"results_{min_f1}.csv", index=False)
            filtered_df.to_excel(output_dir / f"results_{min_f1}.xslx", index=False)
            filtered_df.to_html(output_dir / f"results_{min_f1}.html", index=False)
            filtered_df.describe().to_csv(output_dir / f"summary_{min_f1}.csv", index=False)
            filtered_df.describe().to_html(output_dir / f"summary_{min_f1}.html", index=False)
        df_sorted = df.sort_values("f1", ascending=False)
        df_sorted.to_csv(output_dir / "results_sorted.csv", index=False)
        # show top N chemicals
        N = 20
        df_sorted[:N][["chemical_class", "f1"]].to_csv(output_dir / "top_chemicals.csv", index=False)
        df_sorted[:N][["chemical_class", "f1", "false_positives", "false_negatives"]].to_csv(output_dir / "top_chemicals_failed.csv", index=False)
        df_sorted[:N][["chemical_class", "f1", "false_positives", "sample_false_negatives"]].to_html(output_dir / "top_chemicals_failed.html", index=False)
        # show bottom N chemicals
        df_sorted[-N:][["chemical_class", "f1"]].to_csv(output_dir / "bottom_chemicals.csv", index=False)


if __name__ == "__main__":
    app()