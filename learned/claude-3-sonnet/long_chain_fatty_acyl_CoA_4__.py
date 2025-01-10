"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA core structure - adenine base
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety of CoA"
    
    # Check for thioester linkage (R-C(=O)-S-R)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester linkage"
    
    # Check for phosphate groups with negative charges
    phosphate_pattern = Chem.MolFromSmarts("[OP]([O-])([O-])=O")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 1:
        return False, "Missing phosphate groups with correct charge state"
    
    # Check for CoA structure - pantetheine portion
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine portion of CoA"
        
    # Count carbons in longest chain
    carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            carbons.append(atom.GetIdx())
            
    # Calculate molecular properties
    n_carbons = len([a for a in mol.GetAtoms() if a.GetAtomicNum() == 6])
    n_double_bonds = rdMolDescriptors.CalcNumDoubleBonds(mol)
    
    # Identify features
    features = []
    
    # Check for unsaturation
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        features.append("unsaturated")
        
    # Check for hydroxyl groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2,CH]-[OH]")):
        features.append("hydroxylated")
        
    # Check for keto groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("CC(=O)C")):
        features.append("keto")
        
    # Check for epoxy groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C1OC1")):
        features.append("epoxy")
    
    feature_str = " and ".join(features) if features else "saturated"
    
    # Verify long chain - CoA has about 23 carbons itself, so we need more than that
    # for a long chain fatty acid
    if n_carbons < 35:  # 23 (CoA) + 12 (minimum for long chain)
        return False, f"Carbon count ({n_carbons}) too low for long-chain fatty acyl-CoA"
    
    # Check for correct charge state
    n_negative_charges = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    if n_negative_charges != 4:
        return False, f"Incorrect charge state (expected -4, found {n_negative_charges})"
    
    # Calculate approximate fatty acid chain length
    chain_length = n_carbons - 23  # Subtract CoA carbons
    
    return True, f"Long-chain ({chain_length} carbons) {feature_str} fatty acyl-CoA(4-) with {n_double_bonds} double bonds"