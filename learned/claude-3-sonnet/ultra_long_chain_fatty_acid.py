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

    # Look for a long carbon chain - allow for branching and cyclopropyl groups
    # Match any combination of single/double bonds between carbons
    chain_pattern = Chem.MolFromSmarts('[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]')
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No continuous long carbon chain found"

    # Check molecular weight - should be significant for ultra-long chain
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:  # Rough minimum for C28 fatty acid
        return False, f"Molecular weight ({mol_wt:.1f}) too low for ultra-long chain fatty acid"

    # Count oxygens - should have at least 2 (from COOH) but allow for modifications
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:
        return False, "Too few oxygen atoms for a fatty acid"
    
    # Allow for various modifications (hydroxy, methoxy groups)
    # Check that the majority of atoms are carbons
    total_atoms = len(mol.GetAtoms())
    if carbon_count < total_atoms * 0.6:  # At least 60% should be carbon
        return False, "Carbon content too low for a fatty acid"

    # Check for reasonable molecular formula
    # Allow for various modifications but ensure it's primarily hydrocarbon
    if rdMolDescriptors.CalcMolFormula(mol).count('C') < 28:
        return False, "Insufficient carbon content for ultra-long chain fatty acid"

    return True, f"Ultra-long chain fatty acid with {carbon_count} carbons"