"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid is derived from a diterpene, usually has a C20 skeleton with possible modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 30:  # Allow some flexibility around 20 carbons
        return False, "Carbon count not typical for a diterpenoid"

    # Check for key functional groups common in diterpenoids
    key_groups = ["C(=O)O", "C(=O)", "[OH]", "C=C"]
    for group in key_groups:
        group_pattern = Chem.MolFromSmarts(group)
        if mol.HasSubstructMatch(group_pattern):
            return True, "Contains key functional groups of diterpenoids"
            
    # Consider rearranged C20 skeleton (a more complex analysis may be required)
    # For simplicity, we assume typical diterpenoids have significant hydrocarbon backbone
    # Here, just a generic carbon chain or ring will be assumed
    carbon_skeleton_pattern = Chem.MolFromSmarts("C~C~C~C~C~C~C~C~C~C")
    if mol.HasSubstructMatch(carbon_skeleton_pattern):
        return True, "Matches a pattern consistent with a diterpenoid skeleton"
    
    return False, "Does not meet diterpenoid criteria"

# Metadata for further details
__metadata__ = {
    'chemical_class': {
        'name': 'diterpenoid',
        'definition': 'Any terpenoid derived from a diterpene.',
    }
}