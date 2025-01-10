"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: essential fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are polyunsaturated fatty acids that must be obtained through diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:  # Minimum carbon count for essential fatty acids
        return False, f"Carbon chain too short ({c_count} carbons)"
    if c_count > 40:  # Maximum reasonable length
        return False, f"Carbon chain too long ({c_count} carbons)"

    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    n_double_bonds = len(double_bond_matches)
    
    if n_double_bonds < 2:
        return False, f"Not polyunsaturated (only {n_double_bonds} double bonds)"
    if n_double_bonds > 6:
        return False, f"Too many double bonds ({n_double_bonds})"

    # Check for long continuous carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("[CH2,CH3]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long continuous carbon chain found"

    # Verify it's a fatty acid by checking molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if not ('C' in formula and 'O2' in formula and 'H' in formula):
        return False, "Molecular formula not consistent with fatty acid"

    # Check for phospholipid pattern - if present, extract just the fatty acid part
    phospholipid_pattern = Chem.MolFromSmarts("[P](=[O])([O-])")
    if mol.HasSubstructMatch(phospholipid_pattern):
        return True, "Essential fatty acid (as part of phospholipid)"

    # Calculate the ratio of hydrogens to carbons
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    hc_ratio = h_count / c_count
    
    # Typical H:C ratio for polyunsaturated fatty acids is between 1.5 and 2.0
    if not (1.5 <= hc_ratio <= 2.0):
        return False, f"H:C ratio ({hc_ratio:.1f}) not typical for essential fatty acids"

    return True, f"Polyunsaturated fatty acid with {n_double_bonds} double bonds and {c_count} carbons"