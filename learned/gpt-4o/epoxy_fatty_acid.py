"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a heterocyclic fatty acid containing an epoxide ring as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for epoxide ring pattern, which is a three-membered ring with an oxygen
    # Ensure it forms a ring using a more precise SMARTS
    epoxide_ring_pattern = Chem.MolFromSmarts("C1OC1")
    if not mol.HasSubstructMatch(epoxide_ring_pattern):
        return False, "No epoxide ring found"

    # Check for a long carbon chain characteristic of fatty acids
    # Include unsaturations and allow some branching
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 12 or h_count < (c_count * 2 - 2):  # Allowing for unsaturations typical in fatty acids
        return False, "Too few carbons or incorrect C:H ratio for a fatty acid"

    # Find carboxylic acid group typically at/near chain end
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Ring systems count - should avoid complex molecules with multiple distinct rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 2:
        return False, "Too many rings for a typical epoxy fatty acid"

    # Run more checks if needed depending on misclassifications here noted in debugs
    
    return True, "Contains an epoxide ring and typical characteristics of an epoxy fatty acid (long carbon chain and carboxylic acid group)"