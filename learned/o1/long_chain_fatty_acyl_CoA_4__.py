"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:77989 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    Checks for the presence of a CoA moiety, a thioester linkage, and a long acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple:
            bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
            str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define CoA SMARTS pattern (simplified and corrected)
    coa_smarts = "NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1OC(C(O)C1OP(=O)([O-])[O-])n2cnc3c(ncnc23)N"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Invalid CoA SMARTS pattern"

    # Check for CoA moiety
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Define thioester linkage pattern: C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Find the thioester carbon (start of acyl chain)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"
    thioester_carbon_idx = thioester_matches[0][0]  # Carbonyl carbon

    # Traverse the acyl chain from the carbonyl carbon
    acyl_chain_length = 0
    visited_atoms = set()
    atoms_to_visit = [thioester_carbon_idx]

    while atoms_to_visit:
        atom_idx = atoms_to_visit.pop()
        if atom_idx in visited_atoms:
            continue
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            # Check if the atom is part of the CoA moiety
            if not atom.HasSubstructMatch(coa_pattern):
                acyl_chain_length += 1
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited_atoms:
                        neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
                        # Exclude sulfur atom (S) and atoms in CoA moiety
                        if neighbor_atom.GetAtomicNum() != 16 and not neighbor_atom.HasSubstructMatch(coa_pattern):
                            atoms_to_visit.append(neighbor_idx)

    # Adjust for the carbonyl carbon counted
    acyl_chain_length -= 1  # Exclude the carbonyl carbon from chain length

    if acyl_chain_length < 12:
        return False, f"Acyl chain too short ({acyl_chain_length} carbons), need at least 12"

    return True, f"Contains CoA moiety with long acyl chain of {acyl_chain_length} carbons"