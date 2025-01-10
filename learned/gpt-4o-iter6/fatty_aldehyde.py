"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde with a long aliphatic carbon chain derived from fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for terminal aldehyde group [CH]=O pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1;!R]=O")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No terminal aldehyde group found"
    
    # Look for long aliphatic carbon chain
    chain_length = 0
    for match in aldehyde_matches:
        # Traverse the structure from aldehyde and count consecutive aliphatic carbons
        start_atom_idx = match[0]  # Assuming terminal bond to aldehyde is [C]=O
        visited = set()

        def traverse(atom_idx, length):
            nonlocal chain_length
            if atom_idx in visited:
                return
            visited.add(atom_idx)
            if length > chain_length:
                chain_length = length
            for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    next_idx = neighbor.GetIdx()
                    traverse(next_idx, length + 1)

        traverse(start_atom_idx, 1)
    
    if chain_length < 6:
        return False, f"Carbon chain too short for typical fatty aldehyde (found {chain_length} carbons)"
    
    # Avoid false positives: check for presence of aromatic rings
    if mol.HasSubstructMatch(Chem.MolFromSmarts("a")):
        return False, "Aromatic rings present, not a simple fatty aldehyde"

    return True, "Valid fatty aldehyde: Terminal aldehyde group with a long aliphatic chain"