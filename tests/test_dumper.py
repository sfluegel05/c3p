from c3p.datamodel import Config, EvaluationResult
from c3p.dumper import write_program
from tests.conftest import OUTPUT_DIR


def test_write_program(eval_result_metal_atom: EvaluationResult):
    assert eval_result_metal_atom.train_results.best_result.chemical_class.id == "CHEBI:33521"
    write_program(eval_result_metal_atom, OUTPUT_DIR, Config())
