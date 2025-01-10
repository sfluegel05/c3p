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
    
    # Check for CoA core structure
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety of CoA"
    
    # Check for thioester linkage (R-C(=O)-S-R)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester linkage"
    
    # Check for 4 negative charges (phosphate groups)
    negative_charge_pattern = Chem.MolFromSmarts("[O-]")
    charge_matches = len(mol.GetSubstructMatches(negative_charge_pattern))
    if charge_matches < 4:
        return False, f"Found only {charge_matches} negative charges, need at least 4"
    
    # Check for phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])[O-]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing phosphate groups"
    
    # Calculate carbon chain length in fatty acid portion
    # First, get the total carbon count
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # CoA itself contains 21 carbons, so subtract to get fatty acid chain length
    fatty_acid_carbons = carbon_count - 21
    if fatty_acid_carbons < 12:
        return False, f"Fatty acid chain too short ({fatty_acid_carbons} carbons)"
    
    # Look for long carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long carbon chain found"
    
    # Additional checks for expected CoA features
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine portion of CoA"
    
    ribose_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "Missing ribose portion of CoA"
    
    # If we get here, it's likely a long-chain fatty acyl-CoA(4-)
    features = []
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        features.append("unsaturated")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("CO")):
        features.append("hydroxylated")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("CC(=O)C")):
        features.append("keto")
        
    feature_str = " and ".join(features) if features else "saturated"
    
    return True, f"Long-chain ({fatty_acid_carbons} carbons) {feature_str} fatty acyl-CoA(4-) with thioester linkage"