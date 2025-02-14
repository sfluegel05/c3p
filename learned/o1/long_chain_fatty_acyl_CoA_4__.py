"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:77989 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

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

    # Define CoA SMARTS pattern (simplified)
    coa_smarts = """
    O[C@H]1[C@H](O)[C@H](OP(=O)([O-])OP(=O)([O-])OC[C@H]2O[C@H](n3cnc4c(ncnc43)N)
    [C@H]2OP(=O)([O-])[O-])[C@@H](O)[C@@H]1OP(=O)([O-])[O-]
    """
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
    thioester_atom_idx = thioester_matches[0][0]  # Carbonyl carbon

    # Traverse the acyl chain from the carbonyl carbon
    acyl_chain_length = 0
    visited_atoms = set()
    atoms_to_visit = [thioester_atom_idx]

    while atoms_to_visit:
        atom_idx = atoms_to_visit.pop()
        if atom_idx in visited_atoms:
            continue
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            acyl_chain_length += 1
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited_atoms:
                    # Avoid going back towards the thioester linkage (S atom)
                    if neighbor.GetAtomicNum() != 16:  # Exclude sulfur
                        atoms_to_visit.append(neighbor_idx)

    # Adjust for the carbonyl carbon counted
    acyl_chain_length -= 1  # Exclude the carbonyl carbon from chain length

    if acyl_chain_length < 12:
        return False, f"Acyl chain too short ({acyl_chain_length} carbons), need at least 12"

    return True, f"Contains CoA moiety with long acyl chain of {acyl_chain_length} carbons"