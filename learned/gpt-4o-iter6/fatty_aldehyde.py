"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde with a carbonyl group at one end of an aliphatic carbon chain.

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
    
    # Check for terminal aldehyde group [C]=O pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1]=O")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No terminal aldehyde group found"
    
    # Check if surrounds are aliphatic and chainlike
    for match in aldehyde_matches:
        start_atom_idx = match[0]  # Assuming terminal bond to aldehyde is [C]=O
        visited = set()
        carbon_chain = 0

        def is_aliphatic(atom):
            return atom.GetAtomicNum() == 6 and not atom.IsInRing()

        def traverse(atom_idx):
            nonlocal carbon_chain
            if atom_idx not in visited:
                visited.add(atom_idx)
                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    if is_aliphatic(neighbor):
                        carbon_chain += 1
                        traverse(neighbor.GetIdx())

        traverse(start_atom_idx)

        # Adjust chain length criteria
        if carbon_chain >= 5:
            return True, "Valid fatty aldehyde: Terminal aldehyde group with a suitable aliphatic chain"
    
    return False, f"Carbon chain too short for typical fatty aldehyde (found {carbon_chain} carbons)"