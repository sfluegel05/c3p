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
    
    # Check for phosphate groups - more flexible pattern
    phosphate_patterns = [
        Chem.MolFromSmarts("[OP]([O-])([O-])=O"),
        Chem.MolFromSmarts("[OP](=O)([O-])[O-]"),
        Chem.MolFromSmarts("[P](=O)([O-])([O-])O"),
    ]
    
    has_phosphates = False
    for pattern in phosphate_patterns:
        if mol.HasSubstructMatch(pattern):
            has_phosphates = True
            break
            
    if not has_phosphates:
        return False, "Missing phosphate groups"
    
    # Check for CoA structure - pantetheine portion
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine portion of CoA"
        
    # Calculate molecular properties
    n_carbons = len([a for a in mol.GetAtoms() if a.GetAtomicNum() == 6])
    n_double_bonds = rdMolDescriptors.CalcNumDoubleBonds(mol)
    
    # Check for fatty acid chain
    fatty_chain_pattern = Chem.MolFromSmarts("CCCCCCC")  # At least 7 carbons in chain
    if not mol.HasSubstructMatch(fatty_chain_pattern):
        return False, "Missing long fatty acid chain"
    
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
    
    # Lower carbon threshold - CoA has ~23 carbons, long chain fatty acids typically >10
    if n_carbons < 30:  # Lowered from 35
        return False, f"Carbon count ({n_carbons}) too low for long-chain fatty acyl-CoA"
    
    # More flexible charge check - look for total negative charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge > -3:  # Allow some flexibility in charge representation
        return False, f"Total negative charge ({total_charge}) insufficient for CoA(4-)"
    
    # Calculate approximate fatty acid chain length
    chain_length = n_carbons - 23  # Subtract CoA carbons
    
    return True, f"Long-chain ({chain_length} carbons) {feature_str} fatty acyl-CoA(4-) with {n_double_bonds} double bonds"