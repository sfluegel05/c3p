"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: CHEBI:174450 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde derived from a fatty acid, with a carbonyl group at one end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for exactly one aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if len(aldehyde_matches) != 1:
        return False, f"Found {len(aldehyde_matches)} aldehyde groups (needs exactly 1)"
    aldehyde_carbon = aldehyde_matches[0][0]

    # Verify terminal position (only one adjacent carbon)
    adj_carbons = [atom for atom in mol.GetAtomWithIdx(aldehyde_carbon).GetNeighbors() if atom.GetAtomicNum() == 6]
    if len(adj_carbons) != 1:
        return False, "Aldehyde not terminal"
    first_chain_carbon = adj_carbons[0].GetIdx()

    # Check chain is aliphatic (no aromatic rings)
    if mol.GetAtomWithIdx(first_chain_carbon).GetIsAromatic():
        return False, "Aldehyde attached to aromatic system"

    # Check for forbidden elements (only C, H, O allowed)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Contains non-C/H/O atoms"

    # Check oxygen atoms are only in aldehyde or hydroxyl groups
    bad_oxygen = Chem.MolFromSmarts("[O;!$(O=[CX3H1]);!$(O[H])]")
    if mol.HasSubstructMatch(bad_oxygen):
        return False, "Contains non-aldehyde/hydroxyl oxygen"

    # Check minimum chain length (total carbons >=4)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Chain too short (min 4 carbons)"

    return True, "Terminal aldehyde in aliphatic chain with no forbidden groups"