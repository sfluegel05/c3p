
from pydantic import BaseModel, Field
from typing import List, Optional, Tuple


class ChemicalStructure(BaseModel):
    """Represents a chemical entity with a known specific structure/formula."""
    name: str = Field(..., description="rdfs:label of the structure in CHEBI")
    smiles: str = Field(..., description="SMILES string derived from CHEBI")


class ChemicalClass(BaseModel):
    """Represents a class/grouping of chemical entities."""
    id: str = Field(..., description="id/curie of the CHEBI class")
    name: str = Field(..., description="rdfs:label of the class in CHEBI")
    definition: str = Field(..., description="definition of the structure from CHEBI")
    instances: List[ChemicalStructure] = Field(..., description="positive examples")
    negative_instances: Optional[List[ChemicalStructure]] = Field(default=None, description="negative examples")


class Dataset(BaseModel):
    """
    Represents a dataset of chemical classes.
    """
    classes: List[ChemicalClass]

from typing import Optional


class Config(BaseModel):
    """Experimental setup"""
    llm_model_name: str = "gpt-4o"
    accuracy_threshold: float = 0.5
    max_attempts: int = 3
    max_negative: int = 20


OUTCOME = Tuple[str, Optional[str]]


class Result(BaseModel):
    """Result of running workflow on a chemical class"""
    chemical_class: ChemicalClass
    config: Optional[Config] = None
    code: str
    true_positives: Optional[List[OUTCOME]] = None
    false_positives: Optional[List[OUTCOME]] = None
    true_negatives: Optional[List[OUTCOME]] = None
    false_negatives: Optional[List[OUTCOME]] = None
    attempt: int = 0
    success: bool = True
    best: bool = False
    error: Optional[str] = None
    stdout: Optional[str] = None

    num_true_positives: Optional[int] = None
    num_false_positives: Optional[int] = None
    num_true_negatives: Optional[int] = None
    num_false_negatives: Optional[int] = None

    precision: Optional[float] = None
    recall: Optional[float] = None
    f1: Optional[float] = None

    def calculate(self):
        """Calculate derived statistics"""
        self.num_true_positives = len(self.true_positives or [])
        self.num_false_positives = len(self.false_positives or [])
        self.num_true_negatives = len(self.true_negatives or [])
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
