"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: CHEBI:17525 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate.
    It has a glycerol backbone with a phosphate at position 3 and exactly one acyl group at position 1 or 2.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule fits the criteria
        str: Reason for the decision
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Corrected SMARTS for sn-glycerol-3-phosphate backbone with proper stereochemistry
    # Pattern: [CH2]-O-[C@@H](O)-CH2-O-P(=O)(O)O
    glycerol_pattern = Chem.MolFromSmarts("[CH2]O[C@@H](O)COP(=O)(O)O")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No sn-glycerol-3-phosphate backbone found"
    
    # Get the oxygen atoms attached to C1 and C2 of glycerol from the match
    # Indices in SMARTS: [0:CH2, 1:O (C1-O), 2:C@@H (C2), 3:O (C2-OH), 4:CH2 (C3), 5:O (C3-O-P), 6:P, ...]
    o_c1_idx = matches[0][1]  # O attached to C1
    o_c2_idx = matches[0][3]  # O (hydroxyl) on C2
    
    # SMARTS for ester group: O-C(=O)
    ester_pattern = Chem.MolFromSmarts("[OX2]-C(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check which oxygens (C1-O or C2-O) are part of ester groups
    c1_ester = any(o_c1_idx in match for match in ester_matches)
    c2_ester = any(o_c2_idx in match for match in ester_matches)
    
    ester_count = sum([c1_ester, c2_ester])
    if ester_count != 1:
        return False, f"Found {ester_count} acyl groups, expected 1"
    
    # Check the non-acylated oxygen has a hydroxyl (degree 1 and at least one H)
    non_acylated_o = o_c2_idx if c1_ester else o_c1_idx
    o_atom = mol.GetAtomWithIdx(non_acylated_o)
    if o_atom.GetDegree() != 1 or o_atom.GetTotalNumHs() < 1:
        return False, "Non-acylated position lacks hydroxyl group"
    
    return True, "Monoacyl-sn-glycerol 3-phosphate structure confirmed"