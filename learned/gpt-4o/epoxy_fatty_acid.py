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
    epoxide_ring_pattern = Chem.MolFromSmarts("[C;R]1[O;R][C;R]1")
    if not mol.HasSubstructMatch(epoxide_ring_pattern):
        return False, "No epoxide ring found"

    # Check for a long carbon chain characteristic of fatty acids
    # More flexible definition that allows for unsaturations, branching, and rings (to some extent)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, "Too few carbons for a fatty acid"

    # Improved definition of the carboxylic functional group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Ring systems count - allow for epoxide and minimal additional rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 3:  # Allow up to 3 rings considering ring support around the epoxide
        return False, "Too many rings for a typical epoxy fatty acid"

    return True, "Contains an epoxide ring and typical characteristics of an epoxy fatty acid (long carbon chain and carboxylic acid group)"