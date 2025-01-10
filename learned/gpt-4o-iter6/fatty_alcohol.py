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
        
    # Check for sp3 hybridized carbon bound to hydroxyl groups (O-H)
    oh_pattern = Chem.MolFromSmarts("[CX4;!$(C=[O,N])][OH1]")
    if not mol.HasSubstructMatch(oh_pattern):
        return False, "No aliphatic hydroxyl group found on a carbon chain"

    # Allow complex aliphatic alcohols to some extent including diols, as long as they meet the length
    # Check if the total number of heteroatoms (non-C, non-H) is reasonable
    heteroatom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [6, 1, 8])  # Allowing oxygen
    if heteroatom_count > 1:  # Beyond expected oxygen
        return False, f"Too many heteroatoms for a fatty alcohol, count: {heteroatom_count}"

    # Avoid presence of undesirable functional groups that indicate esters or non-isolated carbonyls
    # Allow ketones only if they are not adjacent to the alcohol
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O.[OH1]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains non-isolated carbonyl suggesting alternative functionality"

    return True, "Contains carbon chain with hydroxyl group(s), consistent with fatty alcohols"