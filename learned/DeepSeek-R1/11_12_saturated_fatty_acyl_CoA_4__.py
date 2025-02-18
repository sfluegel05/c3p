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
    
    # Check for CoA structure components: adenine, ribose, phosphates, pantetheine
    # Using SMARTS patterns for key parts
    # Pantetheine part: SCCNC(=O)CCNC(=O)...
    pantetheine_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine part of CoA"
    
    # Check for phosphate groups (at least two)
    phosphate_pattern = Chem.MolFromSmarts("[O-]P(=O)([O-])O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, "Insufficient phosphate groups"
    
    # Check for adenine (approximate pattern)
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2N")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety"
    
    # Check for thioester group (S-C(=O)-C)
    thioester_pattern = Chem.MolFromSmarts("[S]-[C](=[O])-[CX4]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"
    
    # Get the acyl chain starting from the thioester's carbonyl carbon
    # Assuming the first thioester match is the relevant one
    s_idx, c_o_idx, c_start = thioester_matches[0]
    
    # Function to traverse the acyl chain and check bond between 11-12
    def check_chain(start_atom_idx):
        current_idx = start_atom_idx
        chain = [current_idx]
        for _ in range(11):  # Need to collect up to the 12th carbon
            atom = mol.GetAtomWithIdx(current_idx)
            next_c = None
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in chain:
                    next_c = neighbor.GetIdx()
                    break
            if next_c is None:
                return False  # Chain too short
            chain.append(next_c)
            current_idx = next_c
        # Check bond between 11th and 12th carbons (indices 10 and 11 in chain)
        bond = mol.GetBondBetweenAtoms(chain[10], chain[11])
        if bond and bond.GetBondType() == Chem.BondType.SINGLE:
            return True
        return False
    
    if not check_chain(c_start):
        return False, "11-12 bond is unsaturated or chain too short"
    
    return True, "11,12-saturated fatty acyl-CoA(4-)"