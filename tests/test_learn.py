import pytest
import yaml

from c3p.datamodel import Dataset, Config
from c3p.learn import learn_program_single_iter, learn_program, \
    get_positive_and_negative_train_instances
from tests.test_extractor import DATASET_PATH, test_get_negative_instances


@pytest.mark.integration
def test_learn_single_iter():
    with open(DATASET_PATH, "r") as f:
        dataset = Dataset.model_validate_json(f.read())
    cc = dataset.get_chemical_class_by_name("thienopyrimidine")
    config = Config(llm_model_name="gpt-4o")
    pos, neg = get_positive_and_negative_train_instances(cc, dataset)
    #cc.train_negative = infer_negative_train_examples(dataset, cc)
    result_iter = learn_program_single_iter(cc, pos, neg, config=config)
    assert result_iter is not None
    # get next result - ensuring max one iteration
    result = next(result_iter)
    assert result is not None
    print(yaml.dump(result.model_dump(), sort_keys=False))

@pytest.mark.integration
def test_learn_full():
    with open(DATASET_PATH, "r") as f:
        dataset = Dataset.model_validate_json(f.read())
    cc = dataset.get_chemical_class_by_name("thienopyrimidine")
    config = Config(llm_model_name="gpt-4o")
    pos, neg = get_positive_and_negative_train_instances(cc, dataset)
    #cc.train_negative = infer_negative_train_examples(dataset, cc)
    rset = learn_program(cc, pos, neg, config=config)
    assert rset is not None
    print(yaml.dump(rset.model_dump(), sort_keys=False))