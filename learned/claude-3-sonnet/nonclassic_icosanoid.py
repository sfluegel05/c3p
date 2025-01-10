"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: CHEBI:78512 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    Nonclassic icosanoids are biologically active signalling molecules made by 
    oxygenation of C20 fatty acids, excluding classic icosanoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
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
    if c_count != 20:
        return False, f"Must have exactly 20 carbons, found {c_count}"

    # Count oxygens (should have at least 2 - one from COOH plus at least one more)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for oxygenated fatty acid"

    # Check for double bonds (should have multiple)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_matches < 1:
        return False, "No carbon-carbon double bonds found"

    # Look for patterns that would indicate classic icosanoids
    
    # Prostaglandin pattern (5-membered ring with two side chains)
    prostaglandin_pattern = Chem.MolFromSmarts("[CH2][CH2][CH]1[CH2][CH2][CH]([CH2][CH2][CH2]C(=O)[OH])[CH]1")
    if mol.HasSubstructMatch(prostaglandin_pattern):
        return False, "Contains prostaglandin core structure"

    # Thromboxane pattern (6-membered ring with oxygen)
    thromboxane_pattern = Chem.MolFromSmarts("C1OCCCC1")
    if mol.HasSubstructMatch(thromboxane_pattern):
        return False, "Contains thromboxane-like ring structure"
    
    # Look for common oxygen-containing functional groups
    hydroxy_pattern = Chem.MolFromSmarts("[OH1]")
    epoxy_pattern = Chem.MolFromSmarts("C1OC1")
    
    hydroxy_count = len(mol.GetSubstructMatches(hydroxy_pattern))
    epoxy_count = len(mol.GetSubstructMatches(epoxy_pattern))
    
    if hydroxy_count == 0 and epoxy_count == 0:
        return False, "No hydroxyl or epoxy groups found"

    # Calculate degree of unsaturation
    double_bonds = rdMolDescriptors.CalcNumDoubleBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    total_unsaturation = double_bonds + rings
    
    if total_unsaturation < 2:
        return False, "Insufficient unsaturation for nonclassic icosanoid"

    # If all checks pass, it's likely a nonclassic icosanoid
    return True, f"C20 oxygenated fatty acid with {hydroxy_count} hydroxyl groups, {epoxy_count} epoxy groups, and {double_bond_matches} double bonds"