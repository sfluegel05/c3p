"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: CHEBI:143872 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    Must have: thioester group, CoA backbone with diphosphate groups, adenine moiety, and total charge -4.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for thioester group (S-C(=O)-)
    thioester_pattern = Chem.MolFromSmarts("[SX2]-[CX3](=[OX1])")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester group found"

    # Check for diphosphate group (O-P-O-P-O with two [O-] charges)
    diphosphate_pattern = Chem.MolFromSmarts("[O][P](=[O])([O-])[O][P](=[O])([O-])[O]")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "Diphosphate group not found"

    # Check for additional phosphate with two [O-] charges (adenine-linked phosphate)
    phosphate_pattern = Chem.MolFromSmarts("[O][P]([O-])([O-])=O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing doubly charged phosphate group"

    # Check for adenine moiety (purine base pattern)
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)N")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Adenine moiety not detected"

    # Verify total charge is -4
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != -4:
        return False, f"Total charge {total_charge} ≠ -4"

    return True, "Contains CoA backbone with thioester, diphosphate groups, and correct charge"