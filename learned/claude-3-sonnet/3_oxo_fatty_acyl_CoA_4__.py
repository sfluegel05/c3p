"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA basic structure with more flexible patterns
    # Adenine pattern that accounts for aromaticity
    adenine = Chem.MolFromSmarts("[n;H0]1c[n;H0]c2c1[n;H0]c[n;H0]c2N")
    # More flexible ribose pattern
    ribose = Chem.MolFromSmarts("OC1OC(CO)C(O)C1O")
    
    if not mol.HasSubstructMatch(adenine):
        return False, "Missing adenine structure"
    if not mol.HasSubstructMatch(ribose):
        return False, "Missing ribose structure"

    # Check for phosphate groups with multiple patterns
    phosphate_patterns = [
        "P(=O)([O-])([O-])",  # Doubly charged phosphate
        "P(=O)([O-])(O)",     # Singly charged phosphate
        "P(=O)(O)(O)"         # Neutral phosphate (could be protonated form)
    ]
    
    total_phosphates = 0
    for pattern in phosphate_patterns:
        phos = Chem.MolFromSmarts(pattern)
        total_phosphates += len(mol.GetSubstructMatches(phos))
    
    if total_phosphates < 3:
        return False, f"Expected at least 3 phosphate groups, found {total_phosphates}"

    # Count negative charges
    negative_charges = sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    if abs(negative_charges) < 4:
        return False, f"Expected at least 4 negative charges, found {abs(negative_charges)}"

    # Check for 3-oxo group pattern with thioester
    oxo_thioester = Chem.MolFromSmarts("C(=O)CC(=O)S")
    if not mol.HasSubstructMatch(oxo_thioester):
        return False, "Missing 3-oxo group with thioester linkage"

    # Check for fatty acid chain - more flexible pattern
    fatty_chain = Chem.MolFromSmarts("C(~C~C~C)[C,$(C=C)]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "Missing sufficient fatty acid chain"

    # Check for pantetheine part with more flexible pattern
    pantetheine = Chem.MolFromSmarts("CCNC(=O)CCNC(=O)[CH](O)C(C)(C)")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing pantetheine moiety"

    # Check for complete CoA linkage pattern
    coa_linkage = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[CH](O)C(C)(C)COP([O-,O])")
    if not mol.HasSubstructMatch(coa_linkage):
        return False, "Missing proper CoA linkage pattern"

    # Verify the 3-oxo structure position
    oxo_pattern = Chem.MolFromSmarts("C(=O)CC(=O)SCCNC")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "3-oxo group not properly positioned relative to CoA"

    return True, "Contains 3-oxo fatty acid chain connected to CoA(4-) via thioester bond with proper phosphate groups"