import pytest
import requests
import yaml

from c3p.clients.chebifier import ChebifierClient
from c3p.datamodel import Dataset, ChemicalClass, ChemicalStructure
from c3p.learn import safe_name
from tests.conftest import OUTPUT_DIR
from tests.test_extractor import DATASET_PATH

CHEBIFIER_BASE_URL = "https://chebifier.hastingslab.org/api"

@pytest.fixture
def client():
    return ChebifierClient()

@pytest.mark.parametrize("smiles,expected", [
    ("CSCCCCCCC(NO)C(O)=O", None),
    ("CCCCSCCCCCCC(NO)C(O)=O", None),
    ("O", None),
    #("[Mc]", None),
    ("OO", None),
])
def test_classify(client, smiles, expected):
    print(f"Testing classify with {smiles}")
    results = client.classify(smiles)
    assert results
    for result in results:
        print(yaml.dump(result.model_dump(), sort_keys=False))

def old_test_classify_dataset(client):
    dataset = Dataset(
        classes = [
            ChemicalClass(
                id = "CHEBI:10036",
                name = "wax ester",
                definition = "",
                all_positive_examples = [
                    "C(CCCCCCCC)CCCCC(OCCCCCCCCCCCCCC)=O",
                    "O(CCCCCCCCCCCCCCCCCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\CCCCC",
                ],
            )
        ],
        structures = [
        ],
        validation_examples = [
            "C(CCCCCCCC)CCCCC(OCCCCCCCCCCCCCC)=O",
            "O(CCCCCCCCCCCCCCCCCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\CCCCC",
            "O[C@@H](CCCCCCC(O)=O)/C=C\\C/C=C\\CCCCC",
        ]
    )
    er = client.classify_dataset(dataset)
    print(yaml.dump(er.model_dump(), sort_keys=False))

def test_classify_dataset(client):
    with open(DATASET_PATH, "r") as f:
        dataset = Dataset.model_validate_json(f.read())
    dataset.classes = dataset.classes[:3]
    examples = set()
    classes = []
    for c in dataset.classes:
        if len(c.all_positive_examples) < 100:
            c.all_positive_examples = c.all_positive_examples[:10]
            classes.append(c)
            examples.update(c.all_positive_examples)
        if len(classes) > 3:
            break
    dataset.validation_examples = list(examples)
    print(f"Num examples: {len(examples)}")
    dataset.classes = classes
    print(f"Num classes: {len(classes)}")
    ers = client.classify_dataset(dataset)
    for er in ers:
        print(yaml.dump(er.model_dump(), sort_keys=False))


def test_classify_dataset_full(client):
    with open("../results/2025/benchmark/dataset.json", "r") as f:
        dataset = Dataset.model_validate_json(f.read())
    ers = client.classify_dataset(dataset)
    output_dir = OUTPUT_DIR / "chebifier"
    output_dir.mkdir(exist_ok=True, parents=True)
    for er in ers:
        sn = safe_name(er.test_result.chemical_class.name)
        with open(output_dir / f"{sn}.json", "w") as f:
            f.write(er.model_dump_json(indent=2))





