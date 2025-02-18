"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI: long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    A long-chain fatty acyl-CoA(4-) has a thioester-linked long fatty acid chain, CoA structure,
    and a formal charge of -4 from deprotonated phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Check for thioester group (S-C(=O))
    thioester_pattern = Chem.MolFromSmarts("[SX2]C(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester group (S-C=O) found"

    # Check pantetheine linkage (S-CC-N-C(=O))
    pantetheine_pattern = Chem.MolFromSmarts("[SX2]CCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Pantetheine linkage (S-CC-N-C=O) not found"

    # Check for two phosphate groups with two deprotonated oxygens each
    phosphate_pattern = Chem.MolFromSmarts("[PX4]([O-])([O-])(=O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, f"Found {len(phosphate_matches)} deprotonated phosphate groups, need at least 2"

    # Calculate acyl chain length from thioester
    c_o_match = mol.GetSubstructMatch(thioester_pattern)
    s_idx = c_o_match[0]
    c_idx = c_o_match[1]

    # Find the acyl chain starting from the carbonyl carbon (excluding S side)
    acyl_start = None
    for neighbor in mol.GetAtomWithIdx(c_idx).GetNeighbors():
        if neighbor.GetIdx() != s_idx:
            acyl_start = neighbor
            break
    if not acyl_start:
        return False, "Acyl chain not found"

    # Traverse the acyl chain to count carbons
    visited = {c_idx, s_idx}
    stack = [acyl_start]
    carbon_count = 0

    while stack:
        atom = stack.pop()
        if atom.GetAtomicNum() == 6:  # Carbon atom
            carbon_count += 1
            visited.add(atom.GetIdx())
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited:
                    stack.append(neighbor)

    if carbon_count < 12:
        return False, f"Acyl chain too short ({carbon_count} carbons), needs ≥12"

    # Verify formal charge is -4
    charge = Chem.GetFormalCharge(mol)
    if charge != -4:
        return False, f"Formal charge is {charge}, expected -4"

    return True, "Long-chain fatty acyl-CoA(4-) with all required structural features"