"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: CHEBI:154527 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    Must have: 
    1. Thioester-linked CoA moiety with 3-oxo group (S-C(=O)-CC(=O)-C... pattern)
    2. Adenine-containing CoA structure
    3. Phosphate groups (including deprotonated forms)
    4. Pantetheine chain
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Critical 3-oxo-thioester pattern: S-C(=O)-CC(=O)-[carbon chain]
    # Ensures the 3-oxo is followed by at least one carbon (excluding oxygen-linked)
    thioester_3oxo = Chem.MolFromSmarts('[S]C(=O)CC(=O)[CX4,CX3]')
    if not mol.HasSubstructMatch(thioester_3oxo):
        return False, "Missing 3-oxo-thioester group (S-C(=O)-CC(=O)-C...)"
    
    # CoA structure verification
    # 1. Adenine pattern (n1cnc2ncnc12 allows any substitution)
    adenine_pattern = Chem.MolFromSmarts('n1cnc2ncnc12')
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine ring"
    
    # 2. Phosphate groups (allow deprotonated forms)
    phosphate_pattern = Chem.MolFromSmarts('[O]P(=O)([O-])[O-]')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        # Check for alternative protonation states
        phosphate_pattern_alt = Chem.MolFromSmarts('[O]P(=O)([O])[O]')
        phosphate_alt_matches = mol.GetSubstructMatches(phosphate_pattern_alt)
        if len(phosphate_alt_matches) < 2:
            return False, f"Insufficient phosphate groups (found {len(phosphate_matches)+len(phosphate_alt_matches)})"
    
    # 3. Pantetheine chain (S-C-C-N-C=O)
    pantetheine_pattern = Chem.MolFromSmarts('[S]-C-C-N-C(=O)')
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine chain (SCCNC=O)"
    
    # Additional check: Verify fatty acid chain length (at least 4 carbons in acyl part)
    # Get all matches for the 3-oxo-thioester pattern
    matches = mol.GetSubstructMatches(thioester_3oxo)
    for match in matches:
        s_idx = match[0]
        s_atom = mol.GetAtomWithIdx(s_idx)
        # Follow the thioester bond to carbonyl carbon
        for neighbor in s_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in neighbor.GetBonds()):
                # This is the carbonyl C adjacent to S
                next_c = [n for n in neighbor.GetNeighbors() if n.GetSymbol() == 'C' and n.GetIdx() != s_idx]
                if not next_c:
                    continue
                next_c = next_c[0]
                # Check if this C is part of CC(=O)
                if any(n.GetSymbol() == 'O' and bond.GetBondType() == Chem.BondType.DOUBLE for n in next_c.GetNeighbors() for bond in n.GetBonds()):
                    # Now check the chain after CC(=O)
                    chain_carbon_count = 0
                    visited = set()
                    stack = [next_c]
                    while stack:
                        atom = stack.pop()
                        if atom.GetSymbol() != 'C' or atom.GetIdx() in visited:
                            continue
                        visited.add(atom.GetIdx())
                        chain_carbon_count += 1
                        # Traverse non-ester carbons
                        for nbr in atom.GetNeighbors():
                            if nbr.GetSymbol() == 'C' and nbr.GetIdx() not in visited and not any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in nbr.GetBonds() if bond.GetOtherAtom(nbr).GetSymbol() == 'O'):
                                stack.append(nbr)
                    # Require at least 2 carbons in the chain (3-oxo + 2 = total 5?)
                    if chain_carbon_count >= 2:
                        return True, "Contains 3-oxo-thioester, CoA structure with adenine, phosphates, and sufficient acyl chain"
    
    return False, "Insufficient fatty acid chain length"