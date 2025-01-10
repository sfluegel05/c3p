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
    A secondary amine has a nitrogen atom bonded to two hydrocarbyl groups and one hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Updated SMARTS pattern for more general secondary amine
    # Nitrogen with two carbon connections (could be any hybridization) and one hydrogen
    secondary_amine_smarts = Chem.MolFromSmarts("[NX3;H1][C,c][C,c]") 

    # Check for matches in the molecule
    if mol.HasSubstructMatch(secondary_amine_smarts):
        # Further filtering: check that the nitrogen is not part of an amide, etc.
        amide_smarts = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
        if mol.HasSubstructMatch(amide_smarts):
            return False, "Nitrogen is part of an amide group"

        return True, "Contains secondary amine pattern (nitrogen with two hydrocarbyl groups and one hydrogen)"

    return False, "No secondary amine pattern found"