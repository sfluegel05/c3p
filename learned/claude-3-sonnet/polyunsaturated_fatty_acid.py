"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: CHEBI:36400 polyunsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid is a fatty acid containing more than one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly one carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[C](=O)[O;H,-]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 1:
        return False, "Must have exactly one carboxylic acid group"
    
    # Check for carbon chain length between 16 and 24
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 16 or n_carbons > 24:
        return False, "Carbon chain length outside typical range (16-24)"
    
    # Look for multiple cis double bonds
    cis_double_bond_pattern = Chem.MolFromSmarts("/C=C/")
    cis_double_bond_matches = mol.GetSubstructMatches(cis_double_bond_pattern)
    if len(cis_double_bond_matches) < 2:
        return False, "Fewer than two cis double bonds found"
    
    # Check for typical polyunsaturated fatty acid double bond positions
    double_bond_positions = []
    for match in cis_double_bond_matches:
        bond = mol.GetBondBetweenAtoms(match[0], match[1])
        double_bond_positions.append(bond.GetBeginAtomIdx())
    
    typical_positions = [9, 12, 15, 18, 21]  # Common double bond positions
    if not any(pos in typical_positions for pos in double_bond_positions):
        return False, "Double bond positions atypical for polyunsaturated fatty acids"
    
    # Exclude other lipid classes
    # Check for absence of phosphate groups (phospholipids)
    phosphate_pattern = Chem.MolFromSmarts("[P]")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Likely a phospholipid (phosphate group present)"
    
    # Check for absence of glycerol backbone (triglycerides)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if mol.HasSubstructMatch(glycerol_pattern):
        return False, "Likely a triglyceride (glycerol backbone present)"
    
    return True, "Meets criteria for polyunsaturated fatty acid"