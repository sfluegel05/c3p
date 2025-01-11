import pytest
from oaklib import get_adapter
from oaklib.datamodels.vocabulary import HAS_DBXREF

from c3p.datamodel import Dataset, Config
from c3p.extractor import db_to_dataframe, create_benchmark, validate_dataset
from c3p.learn import run_code, \
    get_positive_and_negative_train_instances
from tests.conftest import OUTPUT_DIR, INPUT_DIR, TEST_PROGRAM_DIR

DATASET_PATH = OUTPUT_DIR / "dataset.json"

#@pytest.mark.integration
def test_create_benchmark():
    db = "sqlite:obo:chebi"
    chebi = get_adapter(db)
    session = chebi.session
    df = db_to_dataframe(session, "CHEBI")
    assert df is not None
    dataset = create_benchmark(df, session)
    assert dataset is not None
    assert len(dataset.classes) > 0
    assert len(dataset.structures) > 0
    with open(DATASET_PATH, "w") as f:
        f.write(dataset.model_dump_json(indent=2))

@pytest.mark.integration
def test_validate():
    with open(DATASET_PATH, "r") as f:
        dataset = Dataset.model_validate_json(f.read())
    bad = []
    for s, sanitized in validate_dataset(dataset):
        print(f"# BAD\n  ORIGINAL: {s.smiles}\n  SANITIZED: {sanitized}\n  NAME: {s.name}")
        bad.append((s, sanitized))
    assert len(bad) == 0


@pytest.mark.integration
def test_get_negative_instances():
    with open(DATASET_PATH, "r") as f:
        dataset = Dataset.model_validate_json(f.read())
    cc = [c for c in dataset.classes if c.name == "diterpene"].pop()
    pos, neg = get_positive_and_negative_train_instances(cc, dataset)
    assert pos
    assert neg
    assert len(neg) > len(pos)
    # assert all(n not in smiles for n in neg_exs)

@pytest.mark.integration
def test_eval():
    with open(DATASET_PATH, "r") as f:
        dataset = Dataset.model_validate_json(f.read())
    cc = [c for c in dataset.classes if c.name == "thienopyrimidine"].pop()
    config = Config(llm_model_name="lbl/gpt-4o", max_negative_to_test=1000)
    result = learn_and_evaluate_for_class(cc, config, dataset)
    assert result is not None
    with open(OUTPUT_DIR / "result.json", "w") as f:
        f.write(result.model_dump_json(indent=2))

@pytest.mark.integration
def test_run_code():
    with open(DATASET_PATH, "r") as f:
        dataset = Dataset.model_validate_json(f.read())
    prog_path = TEST_PROGRAM_DIR / "thioester.py"
    with open(prog_path, "r") as f:
        code_str = f.read()
    smiles_strs = [s.smiles for s in dataset.structures]
    matches = []
    for smiles, is_cls, reason, md in run_code(code_str, "is_thioester", smiles_strs):
        if is_cls:
            matches.append(smiles)
    assert len(matches) > 0





