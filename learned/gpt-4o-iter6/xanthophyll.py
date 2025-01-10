"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is a carotenoid that is oxygenated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General pattern of polyene chain characteristic in carotenoids
    polyene_pattern = Chem.MolFromSmarts("C=C(-C=C)*")  # A more flexible pattern
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "Basic polyene chain characteristic of carotenoids not found"
    
    # Ensure there is a significant polyene backbone
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2)
    if num_double_bonds < 5:
        return False, "Not enough conjugated double bonds for a carotenoid"

    # Check for the presence of oxygen atoms indicating oxidation
    oxy_func_groups = [
        "[OX2H]",  # hydroxy group
        "[C]=O",   # ketone group
        "[R0, R1; $([O]), $([O-][#6])]", # epoxy and other relevant groups
    ]
    
    oxygen_presence = any(
        mol.HasSubstructMatch(Chem.MolFromSmarts(pattern))
        for pattern in oxy_func_groups
    )
    
    if not oxygen_presence:
        return False, "No significant oxygen-based functional groups found"

    return True, "Contains features of an oxygenated carotenoid (xanthophyll)"