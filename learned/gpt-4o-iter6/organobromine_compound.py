"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compounds
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
    
    # Check for presence of bromine
    contains_bromine = any(atom.GetAtomicNum() == 35 for atom in mol.GetAtoms())
    
    if not contains_bromine:
        return False, "No bromine atom found"
    
    # Enhanced querying for carbon-bromine bonds
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 35:  # Look for Br
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Check if it's bonded to a carbon
                    return True, "Contains a carbon-bromine bond"
    
    return False, "No carbon-bromine bond found"