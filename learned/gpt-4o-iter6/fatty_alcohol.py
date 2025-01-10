"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

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
    
    # Count total carbon atoms in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Special condition: Exclude aromatic systems
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Aromatic system present, not a simple fatty alcohol"
      
    # Check for at least one hydroxyl (-OH) group directly attached to a carbon
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")  # Hydroxyl group pattern for aliphatic carbons
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No aliphatic hydroxyl group found on a carbon chain"

    # Ensure no other functional groups skew the classification, e.g., multiple carbonyls, etc.
    # Simple heuristic: Single non-repeating hydroxyl, few branching or other heteroatom chains
    
    # Analyze the hydroxyl groups
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) > 2:
        return False, "More than two hydroxyl groups, potentially not a simple fatty alcohol"

    # Reject complex structures incorrectly identified (e.g., present in multiple distinct chains or alternative functional groups)
    n_heavy_atoms_exceeding_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (6, 8))  # Allow O (for hydroxyl)
    if n_heavy_atoms_exceeding_carbon > 3:
        return False, "Presence of non-alcohol functionality suggesting alternative classification"
    
    # Check carbon count range
    if carbon_count < 3:
        return False, f"Too few carbon atoms ({carbon_count}), need at least 3"
    if carbon_count > 27 and len(hydroxyl_matches) > 1:
        return False, f"Long chain with multiple hydroxyls suggests structure may not be simple fatty alcohol"

    return True, "Contains a carbon chain with the correct hydroxyl group(s), characteristic of fatty alcohols"