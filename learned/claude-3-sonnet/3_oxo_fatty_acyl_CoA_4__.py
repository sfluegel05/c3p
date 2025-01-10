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
    adenine = Chem.MolFromSmarts("n1c(nc2c(N)ncnc12)")
    ribose = Chem.MolFromSmarts("OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)")
    if not (mol.HasSubstructMatch(adenine) and mol.HasSubstructMatch(ribose)):
        return False, "Missing CoA core structure (adenine/ribose)"

    # Check for phosphate groups and negative charges
    # Match any phosphate group with negative charge(s)
    phosphate = Chem.MolFromSmarts("P([O-,OH])([O-,OH])(=O)[O-,OH]")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate))
    if phosphate_matches != 3:
        return False, f"Expected 3 phosphate groups, found {phosphate_matches}"

    # Count total negative charges
    negative_charges = 0
    for atom in mol.GetAtoms():
        negative_charges += abs(min(0, atom.GetFormalCharge()))
    
    if negative_charges != 4:
        return False, f"Expected 4 negative charges, found {negative_charges}"

    # Check for 3-oxo group and thioester linkage
    oxo_thioester = Chem.MolFromSmarts("C(=O)CC(=O)S")
    if not mol.HasSubstructMatch(oxo_thioester):
        return False, "Missing 3-oxo group with thioester linkage"

    # Check for fatty acid chain (at least 4 carbons)
    fatty_chain = Chem.MolFromSmarts("[CH2,CH3][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "Missing sufficient fatty acid chain"

    # Check for pantetheine part
    pantetheine = Chem.MolFromSmarts("CCNC(=O)CCNC(=O)[CH](O)C(C)(C)")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing pantetheine moiety"

    # Count total carbons in fatty acid portion
    # First find the thioester carbon
    thioester = mol.GetSubstructMatch(Chem.MolFromSmarts("C(=O)S"))
    if not thioester:
        return False, "Cannot locate thioester group"
    
    # Verify it's part of a 3-oxo structure
    oxo_carbon = thioester[0]
    neighbors = [n for n in mol.GetAtomWithIdx(oxo_carbon).GetNeighbors() 
                if n.GetAtomicNum() == 6]
    if not any(n.GetTotalNumHs() == 2 for n in neighbors):
        return False, "Missing proper 3-oxo group structure"

    return True, "Contains 3-oxo fatty acid chain connected to CoA(4-) via thioester bond"