"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: CHEBI:33818 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound (contains at least one carbon-bromine bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for carbon-bromine bond
    br_pattern = Chem.MolFromSmarts("[C]-[Br]")
    
    # Check for at least one C-Br bond
    if mol.HasSubstructMatch(br_pattern):
        return True, "Contains at least one carbon-bromine bond"
    else:
        return False, "No carbon-bromine bonds found"