"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol with a carbon chain of 3 to >27 atoms, which may be saturated/unsaturated and branched/unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromaticity: Fatty alcohols must be purely aliphatic
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Aromatic system present, not a simple fatty alcohol"

    # Count total carbon atoms in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, f"Too few carbon atoms ({carbon_count}), need at least 3"
        
    # Check for aliphatic hydroxyl (-OH) groups
    oh_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    if not mol.HasSubstructMatch(oh_pattern):
        return False, "No aliphatic hydroxyl group found on a carbon chain"

    # Allow complex aliphatic alcohols to some extent including diols, as long as they meet the length
    # Check if the total number of heteroatoms (non-C, non-H) is reasonable
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count > 2 and carbon_count < 8:
        return False, f"Too many heteroatoms for small carbon chain size, oxygen count: {oxygen_count}"

    # Avoid alternative functional classifications like non-isolated ester or ketone/aldehydes
    undesirable_functionality = Chem.MolFromSmarts("[C](=O)[!C]")
    if mol.HasSubstructMatch(undesirable_functionality):
        return False, "Contains non-isolated carbonyl suggesting alternative functionality"

    return True, "Contains carbon chain with hydroxyl group(s), consistent with fatty alcohols"