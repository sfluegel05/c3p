"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: CHEBI:76946 ultra-long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    Ultra-long-chain fatty acids have a chain length greater than C27.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[OH]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"
    
    # Count number of carboxylic acid groups - should only be one
    acid_matches = len(mol.GetSubstructMatches(carboxylic_acid))
    if acid_matches > 1:
        return False, f"Found {acid_matches} carboxylic acid groups, should be exactly 1"

    # Count total carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count <= 27:
        return False, f"Carbon chain too short (C{carbon_count}), needs >C27"

    # Look for a continuous carbon chain
    # This SMARTS pattern matches carbon chains (both saturated and unsaturated)
    chain_pattern = Chem.MolFromSmarts('C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C')
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No continuous long carbon chain found"

    # Check molecular weight - should be significant for ultra-long chain
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:  # Rough minimum for C28 fatty acid
        return False, f"Molecular weight ({mol_wt:.1f}) too low for ultra-long chain fatty acid"

    # Additional checks for reasonable fatty acid composition
    # Count oxygens - should have at least 2 (from COOH) but not too many
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:
        return False, "Too few oxygen atoms for a fatty acid"
    if oxygen_count > 10:  # Arbitrary cutoff for reasonable number of oxygens
        return False, "Too many oxygen atoms for a typical fatty acid"

    # Check for reasonable H/C ratio
    hydrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    hc_ratio = hydrogen_count / carbon_count
    if not (1.0 <= hc_ratio <= 2.1):  # Reasonable range for fatty acids
        return False, f"H/C ratio ({hc_ratio:.1f}) outside typical range for fatty acids"

    return True, f"Ultra-long chain fatty acid with {carbon_count} carbons"