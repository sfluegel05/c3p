"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: ultra-long-chain fatty acid
Definition: Any very long-chain fatty acid which has a chain length greater than C27
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_ultra_long_chain, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
        
    # Count total carbons
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 28:
        return False, f"Only {total_carbons} carbons found, need at least 28"
        
    # Look for long carbon chain
    # Allow for both saturated and unsaturated chains
    chain_pattern = Chem.MolFromSmarts('[C,$(C=C)][C,$(C=C)][C,$(C=C)][C,$(C=C)][C,$(C=C)][C,$(C=C)][C,$(C=C)][C,$(C=C)][C,$(C=C)][C,$(C=C)]')
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long carbon chain found"
        
    # Check for reasonable number of oxygens (carboxylic acid + possible hydroxy groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for fatty acid"
    if o_count > 10:
        return False, "Too many oxygen atoms for typical fatty acid"
        
    # Calculate the fraction of carbons that are sp3 hybridized
    sp3_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CX4]')))
    sp2_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CX3]')))
    if sp3_carbons + sp2_carbons == 0:
        return False, "Invalid carbon hybridization"
    
    # Most carbons should be in aliphatic chain
    aliphatic_fraction = (sp3_carbons + sp2_carbons) / total_carbons
    if aliphatic_fraction < 0.8:
        return False, "Structure not primarily aliphatic chain"
        
    # Success case - provide details about the molecule
    details = []
    details.append(f"Contains {total_carbons} carbons")
    if sp2_carbons > 0:
        details.append(f"has {sp2_carbons//2} unsaturated bonds")
    if o_count > 2:
        details.append(f"has {o_count-2} additional oxygen-containing groups")
        
    return True, "Ultra-long chain fatty acid: " + ", ".join(details)