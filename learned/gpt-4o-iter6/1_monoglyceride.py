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

    # Look for the glycerol backbone with ester linkage specifically at position 1
    # Improved SMARTS patterns take into account chirality and correct linkage
    glycerol_pattern = Chem.MolFromSmarts('[C@@H](CO)(OC(=O)[C])O')  # Matches glycerol with ester linkage at position 1

    if mol.HasSubstructMatch(glycerol_pattern):
        return True, "Contains a glycerol backbone with acyl linkage at position 1"

    return False, "Does not have a glycerol backbone with acyl linkage at position 1"