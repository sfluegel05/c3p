"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is a sugar having one or more alcoholic hydroxy groups replaced
    by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern to find sugar backbone (typically an aldose like glucose)
    sugar_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No hexose or typical sugar structure found"
    
    # SMARTS pattern to find amino substitution (O replaced by [NX3,NX4] i.e., NH, NHR, NR2)
    amino_substitution_pattern = Chem.MolFromSmarts("[C@H]1([NX3,NX4])[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O")
    if mol.HasSubstructMatch(amino_substitution_pattern):
        return True, "Contains sugar backbone with one or more OH replaced by amino groups"

    return False, "Sugar backbone found but no amino substitution detected"