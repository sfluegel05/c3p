"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is defined as a monoglyceride in which the acyl substituent is located at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the glycerol backbone with ester linkage at position 1
    # "OC(=O)C[C@@H](O)CO" or "OC(=O)CC(O)CO" are possible substructures
    glycerol_pattern_1 = Chem.MolFromSmarts("OC(=O)C[C@@H](O)CO")
    glycerol_pattern_2 = Chem.MolFromSmarts("OC(=O)CC(O)CO")

    if mol.HasSubstructMatch(glycerol_pattern_1) or mol.HasSubstructMatch(glycerol_pattern_2):
        return True, "Contains a glycerol backbone with acyl linkage at position 1"

    return False, "Does not have a glycerol backbone with acyl linkage at position 1"