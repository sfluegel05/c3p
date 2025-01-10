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

    # Optimized SMARTS pattern to capture typical secondary amines
    # Nitrogen bonded to exactly two carbon atoms and one hydrogen
    secondary_amine_smarts = Chem.MolFromSmarts("[NX3;!$(*=O);!$(*N=[N,O,S])][C,c][C,c]") 

    # Check for matches in the molecule
    if mol.HasSubstructMatch(secondary_amine_smarts):
        # Further filtering to check that the nitrogen is not part of restricted groups
        amide_smarts = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
        nitroso_smarts = Chem.MolFromSmarts("[NX3][N]=[OX1]")
        if mol.HasSubstructMatch(amide_smarts) or mol.HasSubstructMatch(nitroso_smarts):
            return False, "Nitrogen is part of an amide or nitroso group"

        return True, "Contains secondary amine pattern (nitrogen with two hydrocarbyl groups and one hydrogen)"

    return False, "No secondary amine pattern found"