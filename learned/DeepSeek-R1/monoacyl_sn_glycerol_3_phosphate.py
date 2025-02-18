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
    
    # Check for phosphate group attached to glycerol's C3
    # SMARTS for glycerol with phosphate on C3: C(C(COP(=O)(O)O)O)...
    phosphate_pattern = Chem.MolFromSmarts("[CH2][CH]([OH])[CH2]OP(=O)([OH])[OH]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group on glycerol's C3"
    
    # Check for exactly one ester group (acyl) on C1 or C2
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[O;!H0]=[C]([OH0])[OX2]"))
    # Filter esters connected to glycerol's C1 or C2
    glycerol_esters = 0
    for match in ester_matches:
        ester_oxygen = match[1]  # Assuming the oxygen in the ester is the third atom in the SMARTS
        atom = mol.GetAtomWithIdx(ester_oxygen)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() >= 2:  # Check if connected to glycerol's CH2 or CH
                glycerol_esters += 1
    
    if glycerol_esters != 1:
        return False, f"Found {glycerol_esters} acyl groups, expected 1"
    
    # Check remaining hydroxyl groups on glycerol (positions not acylated)
    hydroxyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() > 0:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() >= 2:  # Glycerol's C1/C2
                    hydroxyl_count += 1
    
    if hydroxyl_count < 1:  # At least one OH should remain (since only one acyl group)
        return False, "Missing hydroxyl groups on glycerol"
    
    return True, "Monoacyl-sn-glycerol 3-phosphate structure confirmed"