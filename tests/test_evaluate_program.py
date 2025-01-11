import pytest
import yaml

from c3p.datamodel import Dataset
from c3p.generator import evaluate_program
from c3p.programs import PROGRAM_DIR
from tests.conftest import THIS_DIR, TEST_PROGRAM_DIR


@pytest.mark.parametrize("chemical_class_name", [
    "triglyceride"
])
def test_evaluate_program(chemical_class_name):
    benchmark_dataset = "236-2-5000"
    input_dir = THIS_DIR / ".." / "inputs"
    with open(f"{input_dir}/bench-{benchmark_dataset}.json") as f:
        dataset = Dataset.model_validate_json(f.read())
    chemical_class = None
    positive_instances = None
    negative_instances = None
    for cc in dataset.classes:
        if cc.name == chemical_class_name:
            chemical_class = cc
            positive_instances = [x for x in chemical_class.positive_instances if "*" not in x.smiles]
            negative_instances = chemical_class.negative_instances
            chemical_class.positive_instances = []
            break
    assert chemical_class
    with open(TEST_PROGRAM_DIR / f"{chemical_class.name}.py") as f:
        code = f.read()
    print(f"NUM POSITIVE INSTANCES: {len(positive_instances)}")
    print(f"NEGATIVE INSTANCES: {negative_instances}")
    result = evaluate_program(code, chemical_class, positive_instances, negative_instances)
    print(yaml.dump(result.model_dump()))