"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-) based on its SMILES.
    The molecule must have a CoA structure with a thioester-linked acyl chain where the 11-12 bond is saturated.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for CoA structure components: pantetheine, phosphates, adenine
    pantetheine_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine part of CoA"
    
    phosphate_pattern = Chem.MolFromSmarts("[O-]P(=O)([O-])O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, "Insufficient phosphate groups"
    
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2N")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety"
    
    # Find thioester group (S-C(=O)-C)
    thioester_pattern = Chem.MolFromSmarts("[S]-[C](=[O])-[CX4]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"
    
    # Iterate through matches to find a valid one with exactly 3 atoms
    valid_match = None
    for match in thioester_matches:
        if len(match) == 3:
            valid_match = match
            break
    if not valid_match:
        return False, "Invalid thioester group structure"
    s_idx, c_o_idx, c_start = valid_match
    
    # Traverse the acyl chain to check the 11-12 bond
    def check_saturation(start_idx):
        current_idx = start_idx
        chain = [current_idx]
        # Traverse 11 carbons (positions 1-11), check bond between 11 and 12
        for _ in range(11):
            atom = mol.GetAtomWithIdx(current_idx)
            next_c = None
            # Prefer following the longest chain (non-branching)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in chain:
                    # Check if this is part of the main chain (not a branch)
                    if next_c is None:  # Take the first candidate
                        next_c = neighbor.GetIdx()
                    else:
                        # Branch detected, might need more complex logic
                        pass
            if next_c is None:
                return False  # Chain too short
            chain.append(next_c)
            current_idx = next_c
        # After collecting 12 carbons (indices 0-11), check bond between 11th and 12th
        if len(chain) < 12:
            return False
        bond = mol.GetBondBetweenAtoms(chain[10], chain[11])
        return bond and bond.GetBondType() == Chem.BondType.SINGLE
    
    if not check_saturation(c_start):
        return False, "11-12 bond is unsaturated or chain too short"
    
    # Verify the charge state (4-)
    # CoA(4-) has three phosphate groups (two with -1 each, and one with -2)
    # This is a simplification; exact charge calculation may require more steps
    # Checking for three [O-] in phosphate groups
    negative_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == -1)
    if negative_oxygen < 4:  # Approximate check
        return False, "Insufficient negative charges for CoA(4-)"
    
    return True, "11,12-saturated fatty acyl-CoA(4-)"