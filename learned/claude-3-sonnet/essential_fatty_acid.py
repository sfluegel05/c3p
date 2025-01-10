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
    Common examples include omega-3 and omega-6 fatty acids.

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
    if c_count < 10:  # Allow shorter chains like C10
        return False, f"Carbon chain too short ({c_count} carbons)"
    if c_count > 40:  # Maximum reasonable length
        return False, f"Carbon chain too long ({c_count} carbons)"

    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    n_double_bonds = len(double_bond_matches)
    
    if n_double_bonds < 1:
        return False, f"Not unsaturated (no double bonds)"
    if n_double_bonds > 6:
        return False, f"Too many double bonds ({n_double_bonds})"

    # Look for common essential fatty acid patterns (omega-3, omega-6)
    omega3_pattern = Chem.MolFromSmarts("CC\C=C/C\C=C/C")
    omega6_pattern = Chem.MolFromSmarts("CCCCC\C=C/C\C=C/C")
    
    has_omega3 = mol.HasSubstructMatch(omega3_pattern)
    has_omega6 = mol.HasSubstructMatch(omega6_pattern)
    
    if not (has_omega3 or has_omega6):
        return False, "Does not match omega-3 or omega-6 fatty acid patterns"

    # Check for conjugated double bonds pattern (common in essential fatty acids)
    conjugated_pattern = Chem.MolFromSmarts("C=CC=C")
    has_conjugated = mol.HasSubstructMatch(conjugated_pattern)

    # Check for methylene-interrupted double bonds (common in essential fatty acids)
    methylene_interrupted = Chem.MolFromSmarts("C=CCC=C")
    has_methylene_interrupted = mol.HasSubstructMatch(methylene_interrupted)

    if not (has_conjugated or has_methylene_interrupted):
        return False, "Lacks typical essential fatty acid double bond patterns"

    # Verify it's a fatty acid by checking molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if not ('C' in formula and 'O2' in formula and 'H' in formula):
        return False, "Molecular formula not consistent with fatty acid"

    # Check for phospholipid pattern - if present, still classify as essential fatty acid
    phospholipid_pattern = Chem.MolFromSmarts("[P](=[O])([O-])")
    if mol.HasSubstructMatch(phospholipid_pattern):
        return True, "Essential fatty acid (as part of phospholipid)"

    # Get total hydrogen count (including implicit hydrogens)
    h_count = sum(atom.GetTotalNumHs() + int(atom.GetNumImplicitHs()) 
                 for atom in mol.GetAtoms())
    
    # Calculate molecular properties
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rotatable_bonds < 3:
        return False, "Too rigid for an essential fatty acid"

    return True, f"Essential fatty acid with {n_double_bonds} double bonds and {c_count} carbons"