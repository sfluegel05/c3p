"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: Long-chain fatty acids (C13-C22)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid (C13-C22) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Get all carboxylic acid carbons
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    # Allow for multiple COOH groups but at least one must be terminal
    terminal_cooh = False
    for match in carboxyl_matches:
        carboxyl_carbon = mol.GetAtomWithIdx(match[0])
        if len([n for n in carboxyl_carbon.GetNeighbors() if n.GetAtomicNum() == 6]) == 1:
            terminal_cooh = True
            break
    
    if not terminal_cooh:
        return False, "No terminal carboxylic acid group found"

    # Count total carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 13:
        return False, f"Too few carbons ({num_carbons}) for long-chain fatty acid"
    if num_carbons > 30:  # Allow some extra carbons for modifications
        return False, f"Too many carbons ({num_carbons}) for long-chain fatty acid"
    
    # Look for long aliphatic chain pattern
    # Match carbons connected by single, double, or triple bonds
    chain_pattern = Chem.MolFromSmarts('C(~C~C~C~C~C~C~C~C~C~C~C~C)')
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No continuous carbon chain of sufficient length"
    
    # Count oxygens (fatty acids typically don't have too many oxygens unless modified)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens > 8:  # Allow for some modifications (hydroxy, epoxy, etc.)
        return False, f"Too many oxygen atoms ({num_oxygens}) for typical fatty acid"
    
    # Check for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:  # Typical range for C13-C22 fatty acids with modifications
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for long-chain fatty acids"
    
    # Additional check for reasonable H/C ratio (allowing for some unsaturation)
    num_hydrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    hc_ratio = num_hydrogens / num_carbons
    if not (1.0 <= hc_ratio <= 2.1):  # Allow for unsaturated and slightly modified structures
        return False, f"Unusual hydrogen/carbon ratio ({hc_ratio:.1f})"
    
    # Count rings - fatty acids typically have 0-1 rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings > 2:
        return False, f"Too many rings ({num_rings}) for typical fatty acid"
        
    return True, f"Long-chain fatty acid with {num_carbons} carbons"