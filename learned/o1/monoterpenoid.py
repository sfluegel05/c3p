"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: monoterpenoid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is any terpenoid derived from a monoterpene (C10 skeleton),
    including compounds in which the C10 skeleton has been rearranged or modified
    by the removal of one or more skeletal atoms (generally methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Monoterpenoids typically have around 10 carbons, but may have fewer due to modifications
    if c_count < 7 or c_count > 11:
        return False, f"Number of carbon atoms ({c_count}) not consistent with a monoterpenoid"
    
    # Check for terpenoid functional groups
    # Common functional groups in monoterpenoids include alcohols, ketones, aldehydes, ethers, carboxylic acids
    functional_group_patterns = {
        "alcohol": "[OX2H]",                # Hydroxyl group
        "ketone": "[CX3](=O)[#6]",          # Ketone group
        "aldehyde": "[CX3H1](=O)[#6]",      # Aldehyde group
        "ether": "[OD2]([#6])[#6]",         # Ether group
        "carboxylic_acid": "C(=O)[OH]",     # Carboxylic acid group
        "ester": "C(=O)O[#6]",              # Ester group
        "epoxide": "[C;R][O;R][C;R]",       # Epoxide ring
        "thiol": "[SX2H]",                  # Thiol group
    }

    # Initialize flag for functional group presence
    has_functional_group = False
    for fg_name, smarts in functional_group_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            has_functional_group = True
            break

    if not has_functional_group:
        return False, "No terpenoid functional groups found"
    
    # Optional: Attempt to detect terpenoid skeleton
    # Due to the diversity and complexity of monoterpenoid structures,
    # detecting a common skeleton is challenging
    # For simplicity, we'll rely on carbon count and functional groups

    return True, f"Contains {c_count} carbons and terpenoid functional groups consistent with a monoterpenoid"