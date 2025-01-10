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
    # Parse SMILES and add explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    mol = Chem.AddHs(mol)
    
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
    if num_carbons > 22:
        # Allow slightly more carbons for modifications, but not too many
        if num_carbons > 26:
            return False, f"Too many carbons ({num_carbons}) for long-chain fatty acid"
    
    # Look for carbon chain pattern - more flexible version
    # Match at least 10 carbons in a row (allowing for branching)
    chain_pattern = Chem.MolFromSmarts('C(~C~C~C~C~C~C~C~C~C)')
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No continuous carbon chain of sufficient length"
    
    # Count oxygens (fatty acids typically don't have too many oxygens unless modified)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens > 12:  # Increased limit to allow for more modifications
        return False, f"Too many oxygen atoms ({num_oxygens}) for typical fatty acid"
    
    # Check for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 180 or mol_wt > 600:  # Expanded range to allow for modifications
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for long-chain fatty acids"
    
    # Count rings - fatty acids typically have 0-2 rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings > 2:
        return False, f"Too many rings ({num_rings}) for typical fatty acid"
    
    # Check that most atoms are carbon and hydrogen
    total_atoms = mol.GetNumAtoms(onlyExplicit=False)
    c_and_h = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [1, 6])
    if c_and_h / total_atoms < 0.7:  # At least 70% should be C or H
        return False, "Too many non-carbon/hydrogen atoms for typical fatty acid"
        
    return True, f"Long-chain fatty acid with {num_carbons} carbons"