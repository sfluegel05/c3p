"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:secondary_amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is characterized by having a nitrogen atom bonded to two hydrocarbyl groups and one hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for secondary amine: nitrogen with two carbon groups and one hydrogen
    secondary_amine_smarts = Chem.MolFromSmarts("[NX3;H1;R0][CX4][CX4]") 

    # Check for matches in the molecule
    if mol.HasSubstructMatch(secondary_amine_smarts):
        return True, "Contains secondary amine pattern (nitrogen with two hydrocarbyl groups and one hydrogen)"

    return False, "No secondary amine pattern found"