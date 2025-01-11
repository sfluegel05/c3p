"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: CHEBI:39277 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
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

    # Search for carbon-bromine bonds in any context (aliphatic or aromatic)
    c_br_bond_pattern = Chem.MolFromSmarts("[C,c]-[Br]")
    if mol.HasSubstructMatch(c_br_bond_pattern):
        return True, "Contains at least one carbon-bromine bond"
    else:
        return False, "No carbon-bromine bonds found"