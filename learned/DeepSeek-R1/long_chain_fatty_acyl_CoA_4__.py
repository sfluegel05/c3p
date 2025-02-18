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
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on SMILES.
    Criteria:
    - Contains thioester-linked fatty acid (S-C=O)
    - Has pantetheine linkage (S-CC-N-C=O)
    - Contains two phosphate groups with total 4 deprotonated oxygens (charge -4)
    - Fatty acid chain length ≥12 carbons
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Thioester check (S-C=O)
    thioester = Chem.MolFromSmarts("[SX2]C(=O)")
    if not mol.HasSubstructMatch(thioester):
        return False, "Missing thioester group"

    # Pantetheine linkage (S-CC-N-C=O)
    pantetheine = Chem.MolFromSmarts("[SX2]CCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing pantetheine linkage"

    # Phosphate group validation - must have 2 P atoms with total 4 O- charges
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 2:
        return False, f"Expected 2 P atoms, found {len(p_atoms)}"
    
    total_o_neg = 0
    for p in p_atoms:
        o_neg = sum(1 for n in p.GetNeighbors() 
                   if n.GetAtomicNum() == 8 and n.GetFormalCharge() == -1)
        total_o_neg += o_neg
    
    if total_o_neg != 4:
        return False, f"Found {total_o_neg} deprotonated phosphate oxygens, need 4"

    # Acyl chain length check from thioester
    thio_match = mol.GetSubstructMatch(thioester)
    s_idx = thio_match[0]
    c_idx = thio_match[1]
    
    # Get first carbon in acyl chain (neighbor of carbonyl not S)
    acyl_start = next((n for n in mol.GetAtomWithIdx(c_idx).GetNeighbors() 
                      if n.GetIdx() != s_idx), None)
    if not acyl_start:
        return False, "No acyl chain found"

    # Traverse all connected carbons from start
    visited = set()
    stack = [acyl_start]
    chain_carbons = 0
    
    while stack:
        atom = stack.pop()
        if atom.GetIdx() in visited:
            continue
        visited.add(atom.GetIdx())
        
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_carbons += 1
            # Add non-visited neighbors (exclude backtracking to thioester)
            stack.extend([n for n in atom.GetNeighbors() 
                         if n.GetIdx() not in visited and n.GetIdx() != c_idx])

    if chain_carbons < 12:
        return False, f"Acyl chain too short ({chain_carbons} carbons)"

    # Verify total charge matches -4
    if Chem.GetFormalCharge(mol) != -4:
        return False, "Formal charge not -4"

    return True, "Meets all criteria for long-chain fatty acyl-CoA(4-)"