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
    Criteria: 3-oxo group adjacent to thioester, CoA structure with deprotonated phosphates, total charge -4.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Corrected 3-oxo pattern: R-C(=O)-CH2-C(=O)-S-CoA (beta-keto thioester)
    # Pattern: Two carbonyls separated by CH2, followed by sulfur
    oxo_pattern = Chem.MolFromSmarts("[CX3](=[OX1])-[CX2H2]-[CX3](=[OX1])-[SX2]")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "Missing 3-oxo group adjacent to thioester"
    
    # Check CoA structural components
    # Pantetheine (S-C-C-N-C=O)
    pantetheine = Chem.MolFromSmarts("SCCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing pantetheine in CoA"
    
    # Adenine moiety in adenosine
    adenine = Chem.MolFromSmarts("n1cnc2c([NH2])ncnc12")
    if not mol.HasSubstructMatch(adenine):
        return False, "Missing adenine"
    
    # Verify total charge is -4 (4- charge state)
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != -4:
        return False, f"Total charge {total_charge} != -4"
    
    # Check presence of phosphate groups (at least 3 phosphorus atoms)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 3:
        return False, f"Only {p_count} phosphorus atoms, need ≥3"
    
    return True, "3-oxo group, CoA structure, and charge -4 confirmed"