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
    2. Adenine-containing CoA structure with proper substitution
    3. Phosphate groups (including different protonation states)
    4. Pantetheine chain
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Sanitize to ensure proper aromaticity perception
    Chem.SanitizeMol(mol)
    
    # Critical 3-oxo-thioester pattern: S-C(=O)-CC(=O)-[carbon chain]
    thioester_3oxo = Chem.MolFromSmarts('[SX2]-C(=O)-CC(=O)-[CX4,CX3]')
    if not mol.HasSubstructMatch(thioester_3oxo):
        return False, "Missing 3-oxo-thioester group (S-C(=O)-CC(=O)-C...)"
    
    # CoA structure verification
    # 1. Adenine pattern with NH group (matches actual CoA structure)
    adenine_pattern = Chem.MolFromSmarts('n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine ring with proper substitution"
    
    # 2. Phosphate groups (any protonation state)
    phosphate_pattern = Chem.MolFromSmarts('[PX4](=O)(O)(O)')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need ≥2"
    
    # 3. Pantetheine chain (S-C-C-N-C=O)
    pantetheine_pattern = Chem.MolFromSmarts('[S]-C-C-N-C(=O)')
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine chain (SCCNC=O)"
    
    # Verify fatty acid chain length (minimum 4 carbons in acyl part)
    # Get the thioester carbon chain after CC(=O)
    thioester_match = mol.GetSubstructMatch(thioester_3oxo)
    if not thioester_match:
        return False, "Thioester pattern match failed"
    
    # The CC(=O) is at positions 2 and 3 in the thioester pattern [S]-C(=O)-CC(=O)-...
    # Starting from the last C in CC(=O) (position 3), follow the chain
    chain_start = thioester_match[3]
    chain_carbon_count = 0
    
    # Traverse connected carbons that are part of the fatty acid chain
    visited = set()
    stack = [mol.GetAtomWithIdx(chain_start)]
    while stack:
        atom = stack.pop()
        if atom.GetAtomicNum() != 6 or atom.GetIdx() in visited:
            continue
        visited.add(atom.GetIdx())
        chain_carbon_count += 1
        
        # Add neighboring carbons not part of ester/ketone groups
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and not any(bond.GetBondType() == Chem.BondType.DOUBLE 
                                                 for bond in nbr.GetBonds() 
                                                 if bond.GetOtherAtom(nbr).GetAtomicNum() == 8):
                stack.append(nbr)
    
    if chain_carbon_count < 4:
        return False, f"Fatty acid chain too short ({chain_carbon_count} carbons)"
    
    return True, "Contains 3-oxo-thioester, CoA structure with adenine, phosphates, and sufficient acyl chain"