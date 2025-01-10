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
    
    # Enhanced SMARTS patterns for carbon-bromine bonds
    # "[C]-[Br]" for aliphatic carbon-bromine bonds
    # "[c]-[Br]" for aromatic carbon-bromine bonds
    smarts_patterns = ["[C]-[Br]", "[c]-[Br]"]
    
    # Check each pattern for a match
    for pattern in smarts_patterns:
        cb_bond_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(cb_bond_pattern):
            return True, f"Contains at least one carbon-bromine bond matching pattern: {pattern}"
    
    return False, "No carbon-bromine bond found"