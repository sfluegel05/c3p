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
    if c_count < 16:  # Minimum C16 for essential fatty acids
        return False, f"Carbon chain too short ({c_count} carbons)"
    if c_count > 40:  # Maximum reasonable length
        return False, f"Carbon chain too long ({c_count} carbons)"

    # Count double bonds (using SMARTS that matches any C=C)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if double_bond_pattern is None:
        return None, "Error in SMARTS pattern"
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    n_double_bonds = len(double_bond_matches)
    
    if n_double_bonds < 2:  # Must be polyunsaturated
        return False, f"Not polyunsaturated (needs at least 2 double bonds)"
    if n_double_bonds > 6:
        return False, f"Too many double bonds ({n_double_bonds})"

    # Check for methylene-interrupted double bond system
    # This pattern looks for C=C-C-C=C (double bonds separated by CH2)
    methylene_interrupted = Chem.MolFromSmarts("C=CC[CH2]C=C")
    if methylene_interrupted is None:
        return None, "Error in SMARTS pattern"
    if not mol.HasSubstructMatch(methylene_interrupted):
        return False, "Double bonds not methylene-interrupted"

    # Check for straight carbon chain
    # Most carbons should have 2 connections (middle of chain) or 1 (ends)
    branched = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            if atom.GetDegree() > 3:  # More than 3 connections
                branched = True
                break
    if branched:
        return False, "Contains branched chains"

    # Count non-carbon/hydrogen atoms
    other_atoms = sum(1 for atom in mol.GetAtoms() 
                     if atom.GetAtomicNum() not in [1, 6, 8])
    if other_atoms > 0:
        return False, "Contains atoms other than C, H, O"

    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:  # Typical range for fatty acids
        return False, f"Molecular weight {mol_wt:.1f} outside typical range"

    return True, (f"Polyunsaturated fatty acid with {n_double_bonds} methylene-interrupted "
                 f"double bonds and {c_count} carbons")