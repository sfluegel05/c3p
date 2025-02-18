"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES.
    Criteria: 3-oxo group, CoA structure with deprotonated phosphates, total charge -4.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Check for 3-oxo group (R-C(=O)-C(=O)-S-)
    oxo_pattern = Chem.MolFromSmarts("[CX3]=[OX1]-[CX3]=[OX1]-S")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "Missing 3-oxo group adjacent to thioester"
    
    # Check pantetheine part (S-C-C-N-C=O)
    pantetheine = Chem.MolFromSmarts("SCCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing pantetheine in CoA"
    
    # Check adenine in adenosine
    adenine = Chem.MolFromSmarts("n1cnc2c([NH2])ncnc12")
    if not mol.HasSubstructMatch(adenine):
        return False, "Missing adenine moiety"
    
    # Verify phosphate groups have 4- charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != -4:
        return False, f"Total charge {total_charge} != -4"
    
    # Check deprotonated phosphates (at least three [O-] attached to P)
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O-])[O-]")
    if len(mol.GetSubstructMatches(phosphate_pattern)) < 3:
        return False, "Insufficient deprotonated phosphates"
    
    return True, "3-oxo group, CoA structure, and charge -4 confirmed"