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
    
    # Add hydrogens to accurately count H's on oxygen
    mol = Chem.AddHs(mol)
    
    # SMARTS for sn-glycerol-3-phosphate backbone with correct stereochemistry
    # Matches: [CH2]-O-[C@H](O)-CH2-OP(=O)(O)O
    glycerol_pattern = Chem.MolFromSmarts("[CH2]O[C@H](O)COP(=O)(O)O")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No sn-glycerol-3-phosphate backbone found"
    
    # Get the oxygen atoms attached to C1 and C2 of glycerol
    # Indices in SMARTS pattern: [0:CH2, 1:O (C1-O), 2:C@H (C2), 3:O (C2-OH), 4:C (C3), ...]
    o_c1_idx = matches[0][1]
    o_c2_idx = matches[0][3]
    
    # Function to check if oxygen is part of an ester (O-C=O)
    def is_ester_o(o_atom):
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                for bond in neighbor.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        other_atom = bond.GetOtherAtom(neighbor)
                        if other_atom.GetSymbol() == 'O':
                            return True
        return False
    
    # Check each oxygen for ester linkage
    o_c1 = mol.GetAtomWithIdx(o_c1_idx)
    o_c2 = mol.GetAtomWithIdx(o_c2_idx)
    
    ester_count = 0
    if is_ester_o(o_c1):
        ester_count += 1
    if is_ester_o(o_c2):
        ester_count += 1
    
    if ester_count != 1:
        return False, f"Found {ester_count} acyl groups, expected 1"
    
    # Check the other oxygen has a hydroxyl group
    other_o = o_c2 if is_ester_o(o_c1) else o_c1
    if other_o.GetTotalNumHs() < 1:
        return False, "Non-acylated position lacks hydroxyl group"
    
    return True, "Monoacyl-sn-glycerol 3-phosphate structure confirmed"