"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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

    # Define a simple sugar pattern (a flexible version, considering typical hexoses)
    sugar_pattern = Chem.MolFromSmarts("C[C@@H]([OH])[C@@H]([OH])C(=O)[C@H]([OH])[C@H](O)[C@@H]1[OH]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No typical sugar backbone found"
        
    # Define a pattern for amino substitution, with allowance for some variation
    amino_pattern = Chem.MolFromSmarts("[NX3,NX4][$([C])[$(C(O));!O]]")  # NH-[C] where an OH might be replaced
    if mol.HasSubstructMatch(amino_pattern):
        return True, "Contains sugar backbone with hydroxyl groups replaced by amino groups"

    return False, "Sugar backbone found but no amino substitution detected"