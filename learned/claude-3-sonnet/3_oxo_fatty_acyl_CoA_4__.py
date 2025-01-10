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

    # Check for CoA basic structure
    adenine = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    ribose = Chem.MolFromSmarts("O[CH]1O[CH]([CH](O)[CH]1O)")
    if not (mol.HasSubstructMatch(adenine) and mol.HasSubstructMatch(ribose)):
        return False, "Missing CoA core structure (adenine/ribose)"

    # Check for phosphate groups - using a more general pattern
    # Match phosphate groups with negative charges
    phosphate = Chem.MolFromSmarts("P(=O)([O-])([O-])")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate))
    
    # Also check for alternative phosphate representations
    alt_phosphate = Chem.MolFromSmarts("P(=O)(O)(O)")
    phosphate_matches += len(mol.GetSubstructMatches(alt_phosphate))
    
    if phosphate_matches < 3:
        return False, f"Expected at least 3 phosphate groups, found {phosphate_matches}"

    # Count negative charges more carefully
    negative_charges = 0
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() < 0:
            negative_charges += abs(atom.GetFormalCharge())
    
    if negative_charges < 4:
        return False, f"Expected at least 4 negative charges, found {negative_charges}"

    # Check for 3-oxo group and thioester linkage
    # More specific pattern for 3-oxo-acyl group
    oxo_thioester = Chem.MolFromSmarts("C(=O)CC(=O)S")
    if not mol.HasSubstructMatch(oxo_thioester):
        return False, "Missing 3-oxo group with thioester linkage"

    # Check for fatty acid chain - allow for both saturated and unsaturated chains
    fatty_chain = Chem.MolFromSmarts("[C][C][C][C,$(C=C)]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "Missing sufficient fatty acid chain"

    # Check for pantetheine part
    pantetheine = Chem.MolFromSmarts("CCNC(=O)CCNC(=O)C(O)C(C)(C)")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing pantetheine moiety"

    # Verify the complete structure
    # Count carbons in the fatty acid portion
    thioester = mol.GetSubstructMatch(Chem.MolFromSmarts("C(=O)S"))
    if not thioester:
        return False, "Cannot locate thioester group"
    
    # Verify 3-oxo structure
    oxo_carbon = thioester[0]
    neighbors = [n for n in mol.GetAtomWithIdx(oxo_carbon).GetNeighbors() 
                if n.GetAtomicNum() == 6]
    if not any(n.GetTotalNumHs() == 2 for n in neighbors):
        return False, "Missing proper 3-oxo group structure"

    # Additional check for CoA linkage pattern
    coa_linkage = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])")
    if not mol.HasSubstructMatch(coa_linkage):
        return False, "Missing proper CoA linkage pattern"

    return True, "Contains 3-oxo fatty acid chain connected to CoA(4-) via thioester bond with proper phosphate groups"