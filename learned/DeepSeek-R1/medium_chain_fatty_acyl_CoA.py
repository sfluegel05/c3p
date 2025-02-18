"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:61904 medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA has a CoA moiety linked via a thioester bond to a fatty acid with 6-12 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for thioester group (C(=O)S-C) connected to CoA's sulfur
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2][CX4]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"

    # Verify CoA structure presence (simplified check for pantetheine-phosphate-ribose-adenine)
    # Look for adenine (N-containing aromatic), ribose (O-rich chain), and phosphate groups
    # Simplified check: presence of adenine-like ring and multiple phosphate groups
    adenine_pattern = Chem.MolFromSmarts("n1cnc2ncnc12")
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2-])[OX2]")
    if not mol.HasSubstructMatch(adenine_pattern) or len(mol.GetSubstructMatches(phosphate_pattern)) < 2:
        return False, "Missing CoA components (adenine/phosphate)"

    # Identify the fatty acid chain attached to the thioester
    # Start from the sulfur atom in the thioester and traverse the carbon chain
    try:
        sulfur_idx = thioester_matches[0][2]
        neighbor = [a.GetIdx() for a in mol.GetAtomWithIdx(sulfur_idx).GetNeighbors() if a.GetSymbol() == "C"][0]
    except IndexError:
        return False, "Invalid thioester connectivity"

    # Traverse the carbon chain from the sulfur-connected carbon
    chain = []
    visited = set()
    def traverse(atom_idx, prev_idx):
        if atom_idx in visited:
            return
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == "C":
            chain.append(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() != prev_idx and neighbor.GetSymbol() in ["C", "O", "S"]:
                    traverse(neighbor.GetIdx(), atom_idx)
    traverse(neighbor, sulfur_idx)

    # Calculate chain length (excluding the thioester carbonyl)
    chain_length = len(chain) - 1  # subtract the sulfur-attached carbon
    if not (6 <= chain_length <= 12):
        return False, f"Chain length {chain_length} not in 6-12 range"

    return True, f"Medium-chain ({chain_length} carbons) fatty acyl-CoA"