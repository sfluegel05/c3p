"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    Long-chain fatty acyl-CoAs have:
    - CoA moiety
    - Thioester linkage
    - Fatty acid chain length C13-C22
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define and validate all SMARTS patterns
    patterns = {
        'adenine': 'n1cnc2c(N)ncnc12',
        'thioester': 'C(=O)SC',
        'pantetheine': 'CCNC(=O)CCNC(=O)',
        'diphosphate': 'OP(=O)(O)OP(=O)(O)',
        'ribose_phosphate': 'OCC1OC(n)C(O)C1OP(O)(O)=O'
    }
    
    smarts_patterns = {}
    for name, pattern in patterns.items():
        smarts_mol = Chem.MolFromSmarts(pattern)
        if smarts_mol is None:
            return False, f"Internal error: Invalid SMARTS pattern for {name}"
        smarts_patterns[name] = smarts_mol

    # Check for required CoA structural features
    required_features = ['adenine', 'thioester', 'pantetheine', 'diphosphate']
    for feature in required_features:
        if not mol.HasSubstructMatch(smarts_patterns[feature]):
            return False, f"Missing {feature} moiety required for CoA structure"

    # Find the thioester carbon and traverse the fatty acid chain
    thioester_matches = mol.GetSubstructMatches(smarts_patterns['thioester'])
    if not thioester_matches:
        return False, "Could not identify thioester linkage"
    
    # Get the total number of carbons and oxygens
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # CoA typically has 21-23 carbons depending on configuration
    # Subtract from total to get fatty acid chain length
    # Allow for some variation in CoA structure
    estimated_coa_carbons = 22  # Average
    chain_carbons = total_carbons - estimated_coa_carbons
    
    # Check chain length (C13-C22)
    if chain_carbons < 13:
        return False, f"Fatty acid chain too short (approximately C{chain_carbons}, need C13-C22)"
    if chain_carbons > 22:
        return False, f"Fatty acid chain too long (approximately C{chain_carbons}, need C13-C22)"

    # Additional checks for common features
    # Count unsaturations
    double_bonds = rdMolDescriptors.CalcNumAliphaticDoubleBonds(mol)
    triple_bonds = rdMolDescriptors.CalcNumAliphaticTripleBonds(mol)
    
    # Look for common modifications
    hydroxy_pattern = Chem.MolFromSmarts('CO')
    oxo_pattern = Chem.MolFromSmarts('CC(=O)C')
    
    features = []
    if double_bonds > 0:
        features.append(f"{double_bonds} double bond(s)")
    if mol.HasSubstructMatch(hydroxy_pattern):
        features.append("hydroxy group(s)")
    if mol.HasSubstructMatch(oxo_pattern):
        features.append("oxo group(s)")
    
    feature_str = ", ".join(features)
    if feature_str:
        feature_str = f" with {feature_str}"
    
    return True, f"Long-chain fatty acyl-CoA (approximately C{chain_carbons}){feature_str}"