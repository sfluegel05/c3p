"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA moiety
    # Look for adenine base
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety of CoA"
    
    # Look for ribose phosphate
    ribose_phosphate = Chem.MolFromSmarts("OCC1OC(n2cnc3c2ncnc3)C(O)C1OP(O)(O)=O")
    if not mol.HasSubstructMatch(ribose_phosphate):
        return False, "Missing ribose phosphate portion of CoA"
    
    # Look for pantetheine portion with thioester
    # [CX4]C(=O)NCCC(=O)NCCS represents the pantetheine chain
    # The C(=O)S represents the thioester linkage
    pantetheine_thioester = Chem.MolFromSmarts("[CX4]C(=O)NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(pantetheine_thioester):
        return False, "Missing pantetheine-thioester portion"
    
    # Check for 3-oxo group pattern
    # Pattern: carbon chain - C(=O)CC(=O)S- 
    # This represents the 3-oxo group followed by the thioester
    oxo_pattern = Chem.MolFromSmarts("C-C(=O)CC(=O)S")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "Missing 3-oxo group pattern"
    
    # Additional checks for fatty acid portion
    # Count carbons in the main chain
    # We'll be lenient here as fatty acids can vary in length
    carbon_chain = Chem.MolFromSmarts("CCCCC")  # At least 5 carbons
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "Carbon chain too short for fatty acid portion"

    # Check for reasonable molecular weight
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 800:  # CoA itself is quite large
        return False, "Molecular weight too low for 3-oxo-fatty acyl-CoA"
    
    return True, "Contains CoA moiety, 3-oxo group, and appropriate fatty acid chain"