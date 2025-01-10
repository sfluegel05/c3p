"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compounds
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine based on its SMILES string.
    An organobromine compound contains at least one carbon-bromine bond.

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
    cb_bond_pattern = Chem.MolFromSmarts("[C]-[Br]")
    if mol.HasSubstructMatch(cb_bond_pattern):
        return True, "Contains at least one carbon-bromine bond"
    
    return False, "No carbon-bromine bond found"