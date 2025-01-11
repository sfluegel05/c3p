from copy import copy

from pydantic import BaseModel, Field
from typing import List, Optional, Tuple, Union, Set, Dict

SMILES_STRING = str


# make this hashable?
class ChemicalStructure(BaseModel):
    """Represents a chemical entity with a known specific structure/formula."""
    name: str = Field(..., description="rdfs:label of the structure in CHEBI")
    smiles: SMILES_STRING = Field(..., description="SMILES string derived from CHEBI")

    def __hash__(self):
        return hash(self.smiles)

    def __eq__(self, other):
        if not isinstance(other, ChemicalStructure):
            return NotImplemented
        return self.smiles == other.smiles


class ChemicalClass(BaseModel):
    """Represents a class/grouping of chemical entities."""
    id: str = Field(..., description="id/curie of the CHEBI class")
    name: str = Field(..., description="rdfs:label of the class in CHEBI")
    definition: str = Field(..., description="definition of the structure from CHEBI")
    parents: Optional[List[str]] = Field(default=None, description="parent classes")
    xrefs: Optional[List[str]] = Field(default=None, description="mappings")
    all_positive_examples: List[SMILES_STRING] = []
    #train_positive: List[SMILES_STRING] = []
    #train_negative: Optional[List[SMILES_STRING]] = None  # can be inferred
    #validate_positive: List[SMILES_STRING] = []
    #validate_negative: List[SMILES_STRING] = []
    #num_train_positive: Optional[int] = None
    #num_train_negative: Optional[int] = None
    #num_validate_positive: Optional[int] = None
    #num_validate_negative: Optional[int] = None
    #positive_instances: List[ChemicalStructure] = Field(..., description="positive examples")
    #negative_instances: Optional[List[ChemicalStructure]] = Field(default=None, description="negative examples")

    def lite_copy(self) -> "ChemicalClass":
        """
        Create a copy of the chemical class without the instance fields
        Returns:
        """
        cc = copy(self)
        cc.all_positive_examples = []
        #cc.train_positive = []
        #cc.train_negative = []
        #cc.validate_positive = []
        #cc.validate_negative = []
        return cc

class Dataset(BaseModel):
    """
    Represents a dataset of chemical classes.
    """
    ontology_version: Optional[str] = None
    min_members: Optional[int] = None
    max_members: Optional[int] = None
    classes: List[ChemicalClass]
    structures: List[ChemicalStructure] = None
    validation_examples: Optional[List[SMILES_STRING]] = None

    @property
    def name(self):
        return f"bench-{self.ontology_version}-{self.min_members}-{self.max_members}"

    def all_smiles(self) -> Set[SMILES_STRING]:
        return {s.smiles for s in self.structures}

    def smiles_to_instance(self) -> Dict[SMILES_STRING, ChemicalStructure]:
        return {s.smiles: s for s in self.structures}

    def get_chemical_class_by_id(self, class_id: str) -> ChemicalClass:
        for cc in self.classes:
            if cc.id == class_id:
                return cc
        raise ValueError(f"Class {class_id} not found in dataset")

    def get_chemical_class_by_name(self, class_name: str) -> ChemicalClass:
        for cc in self.classes:
            if cc.name == class_name:
                return cc
        raise ValueError(f"Class {class_name} not found in dataset")

from typing import Optional


class Config(BaseModel):
    """Experimental setup"""
    llm_model_name: str = "gpt-4o"
    experiment_local_name: Optional[str] = None
    f1_threshold: float = 0.8
    max_attempts: int = 4
    min_positive_examples_for_training: int = 23
    min_negative_examples_for_training: int = 23
    min_positive_examples_for_validation: int = 5
    min_negative_examples_for_validation: int = 5
    max_positive_instances: Optional[int] = None
    max_positive_to_test: Optional[int] = None
    max_negative_to_test: Optional[int] = None
    max_positive_in_prompt: int = 50
    max_negative_in_prompt: int = 20
    max_examples_in_feedback: Optional[int] = 25
    test_proportion: float = 0.2
    use_definitions: bool = True

    @property
    def experiment_name(self):
        ln = self.experiment_local_name or "undef"
        model_name = self.llm_model_name
        if "/" in model_name:
            model_name = model_name.split("/")[-1]
        return f"{model_name}-{ln}"


OUTCOME = Union[Tuple[str, Optional[str]], Tuple[str, str, Optional[str]]]

class Outcome(BaseModel):
    smiles: SMILES_STRING
    name: Optional[str] = None
    reason: Optional[str] = None

    def __str__(self):
        return f"{self.smiles} ({self.name}): {self.reason}"

    def __repr__(self):
        return f"{self.smiles} ({self.name}): {self.reason}"


class Result(BaseModel):
    """Result of running workflow on a chemical class"""
    chemical_class: ChemicalClass
    config: Optional[Config] = None
    code: str
    message: Optional[str] = None
    true_positives: Optional[List[Outcome]] = None
    false_positives: Optional[List[Outcome]] = None
    true_negatives: Optional[List[Outcome]] = None
    false_negatives: Optional[List[Outcome]] = None
    sample_true_negatives: Optional[List[Outcome]] = None
    sample_false_negatives: Optional[List[Outcome]] = None
    attempt: int = 0
    reasoning: Optional[str] = None
    success: bool = True  ## True if no runtime errors or compilation errors
    best: bool = False
    error: Optional[str] = None
    stdout: Optional[str] = None

    num_true_positives: Optional[int] = None
    num_false_positives: Optional[int] = None
    num_true_negatives: Optional[int] = None
    num_false_negatives: Optional[int] = None
    num_negatives: Optional[int] = None

    precision: Optional[float] = None
    recall: Optional[float] = None
    f1: Optional[float] = None
    accuracy: Optional[float] = None
    negative_predictive_value: Optional[float] = None

    def calculate(self):
        """Calculate derived statistics"""
        self.num_true_positives = len(self.true_positives or [])
        self.num_false_positives = len(self.false_positives or [])
        if self.num_true_negatives is None:
            self.num_true_negatives = len(self.true_negatives or [])
        if self.num_false_negatives is None:
            self.num_false_negatives = len(self.false_negatives or [])
        if self.num_true_positives + self.num_false_positives:
            self.precision = self.num_true_positives / (self.num_true_positives + self.num_false_positives)
        else:
            self.precision = 0.0
        if self.num_true_positives + self.num_false_negatives:
            self.recall = self.num_true_positives / (self.num_true_positives + self.num_false_negatives)
        else:
            self.recall = 0
        if self.precision and self.recall:
            self.f1 = 2 * (self.precision * self.recall) / (self.precision + self.recall)
        else:
            self.f1 = 0
        if self.num_true_positives + self.num_true_negatives + self.num_false_positives + self.num_false_negatives:
            self.accuracy = (self.num_true_positives + self.num_true_negatives) / (
                    self.num_true_positives + self.num_true_negatives + self.num_false_positives + self.num_false_negatives)
        if self.num_true_negatives + self.num_false_negatives:
            self.negative_predictive_value = self.num_true_negatives / (self.num_true_negatives + self.num_false_negatives)
        else:
            self.negative_predictive_value = 0

class ResultSet(BaseModel):
    """A set of results"""
    best_result: Optional[Result] = None
    results: List[Result]
    sorted_attempts: List[int] = []
    experiment_name: Optional[str] = None

    @classmethod
    def from_results(cls, results: List[Result]) -> "ResultSet":
        """Populate the result set from a list of results"""
        obj = cls(results=results, best_result=max(results, key=lambda r: r.f1))
        if obj.best_result:
            obj.best_result.best = True
        return obj

class EvaluationResult(BaseModel):
    """Result of evaluating a model"""
    train_results: ResultSet
    test_result: Result

    def calculate_reward(self):
        """Calculate derived statistics"""
        tr = self.test_result
        tr.calculate()
        br = self.train_results.best_result
        br.calculate()
        tp_reward = tr.num_true_positives * br.precision
        fp_penalty = tr.num_false_positives * br.precision
        fn_penalty = tr.num_false_negatives * br.negative_predictive_value
        reward = (tp_reward - fp_penalty) - fn_penalty
        return reward


class EvaluationExperiment(BaseModel):
    """Represents an evaluation experiment"""
    config: Config
    evaluation_results: List[EvaluationResult]

class ClassificationResult(BaseModel):
    input_smiles: SMILES_STRING
    class_id: str
    class_name: Optional[str] = None
    is_match: bool
    reason: Optional[str] = None
    confidence: Optional[float] = None
