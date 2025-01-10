"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: CHEBI:37600 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is a monosaccharide with a phosphate group esterified to a hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphate group connected to carbon
    phosphate_patterns = [
        "[CX4]-[OX2]-[P](=[O])([O,OH])([O,OH])",  # Standard phosphate ester
        "[CX4]-[OX2]-P(=O)([O-])([O-])",  # Ionized form
        "[CX4]-[OX2]-P([O-])([O-])=O",    # Alternative ionized form
        "[CX4]-[OX2]-P([OH])(=O)[OH]"     # Fully protonated form
    ]
    
    has_phosphate = False
    for pattern in phosphate_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_phosphate = True
            break
            
    if not has_phosphate:
        return False, "No phosphate ester found"

    # Sugar patterns
    sugar_patterns = [
        # Furanose (5-membered ring)
        "[CH2X4,CH1X4]-1-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[OX2]-1",
        # Pyranose (6-membered ring)
        "[CH2X4,CH1X4]-1-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[OX2]-1",
        # Ribose specific pattern
        "[CH2X4]-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[OX2]",
        # Ketose pattern (fructose-like)
        "[CH2X4](-[OX2])-[CX4](-[OX2])-[CH1X4](-[OX2])-[CH1X4](-[OX2])",
        # Open chain aldose
        "[CH1X4](=O)-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[CH1X4](-[OX2])",
        # Deoxyribose pattern
        "[CH2X4]-[CH1X4]-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[OX2]"
    ]

    has_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_sugar = True
            break

    if not has_sugar:
        return False, "No sugar moiety found"

    # Count hydroxyl/phosphoryl attachment points
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4]-[OX2]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_count < 3:  # Need at least 3 oxygens attached to carbons
        return False, "Insufficient oxygen attachments for a sugar"

    # Verify the structure has appropriate molecular weight
    mol_weight = sum([atom.GetMass() for atom in mol.GetAtoms()])
    if mol_weight < 150:  # Minimum weight for a phosphorylated sugar
        return False, "Molecule too small to be a phospho sugar"

    return True, "Contains monosaccharide with phosphate ester"