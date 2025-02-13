from pathlib import Path

from c3p.benchmark_validator import write_runner_cases_from_eval_dir


def test_create_cases():
    working_dir = Path("../results/2025/ensembl-5/eval/")
    write_runner_cases_from_eval_dir(working_dir)