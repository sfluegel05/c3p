"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a heterocyclic fatty acid containing an epoxide ring.

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
    
    # Look for epoxide ring pattern, a three-membered ring with an oxygen
    epoxide_pattern = Chem.MolFromSmarts("[C]-[O]-[C]")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"

    # Check for a long carbon chain with carboxylic acid group
    # Let's define a fatty acid signature here: long carbon chain with -COOH end group
    # Carbon chain length threshold is arbitrarily set to >12 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, "Too few carbons to be considered a fatty acid"

    # Find carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    return True, "Contains an epoxide ring and characteristics of a fatty acid (long carbon chain and carboxylic acid group)"