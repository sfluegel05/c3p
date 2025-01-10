"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is a carotenoid that is oxygenated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a flexible long conjugated chain pattern
    chain_pattern = Chem.MolFromSmarts("C=C" + "(-,-,[#6])" * 5 + "=C=C")  
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No proper long conjugated chain pattern, which is typical for carotenoids"
    
    # Check for the presence of oxygen atoms in various functional groups
    o_pattern = Chem.MolFromSmarts("[#8]")  # General oxygen pattern to match -OH, =O, ethers, etc.
    if not mol.HasSubstructMatch(o_pattern):
        return False, "Carotene derivative must be oxygenated"

    # Check for presence of at least one six-membered carbon ring
    ring_info = mol.GetRingInfo()
    six_membered_carbon_rings = [
        ring for ring in ring_info.AtomRings() 
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)
    ]
    if len(six_membered_carbon_rings) < 1:
        return False, "Carotenoid backbone typically includes cyclic rings"

    return True, "Contains features of an oxygenated carotenoid (xanthophyll)"