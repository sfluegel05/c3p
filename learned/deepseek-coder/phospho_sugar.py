"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: CHEBI:18179 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

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

    # More flexible phosphate group pattern
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        # Check for alternative phosphate patterns
        phosphate_pattern2 = Chem.MolFromSmarts("[PX4]([OX2])([OX2])[OX2]")
        phosphate_matches2 = mol.GetSubstructMatches(phosphate_pattern2)
        if len(phosphate_matches2) == 0:
            return False, "No phosphate group found"

    # More flexible sugar backbone pattern
    # Matches both linear and cyclic sugars with at least 3 carbons and 2 hydroxyl groups
    sugar_pattern = Chem.MolFromSmarts("[C;H0,H1,H2][C;H0,H1,H2]([OH])[C;H0,H1,H2]([OH])")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        # Check for ring structures with at least two hydroxyl groups
        ring_sugar_pattern = Chem.MolFromSmarts("[C;H0,H1,H2]1[C;H0,H1,H2][C;H0,H1,H2]([OH])[C;H0,H1,H2]([OH])[C;H0,H1,H2]1")
        ring_sugar_matches = mol.GetSubstructMatches(ring_sugar_pattern)
        if len(ring_sugar_matches) == 0:
            # Check for more complex sugar structures
            complex_sugar_pattern = Chem.MolFromSmarts("[C;H0,H1,H2][C;H0,H1,H2]([OH])[C;H0,H1,H2]([OH])[C;H0,H1,H2]")
            complex_sugar_matches = mol.GetSubstructMatches(complex_sugar_pattern)
            if len(complex_sugar_matches) == 0:
                return False, "No sugar backbone found"

    # More comprehensive ester bond pattern
    ester_pattern = Chem.MolFromSmarts("[C;H0,H1,H2][OX2][PX4](=[OX1])([OX2])[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        # Check for alternative ester patterns
        ester_pattern2 = Chem.MolFromSmarts("[C;H0,H1,H2][OX2][PX4]([OX2])([OX2])[OX2]")
        ester_matches2 = mol.GetSubstructMatches(ester_pattern2)
        if len(ester_matches2) == 0:
            return False, "Phosphate not esterified to sugar"

    # Relaxed carbon count restriction to 3-10 to accommodate larger sugars
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3 or c_count > 10:
        return False, "Not a monosaccharide (incorrect number of carbons)"

    # Check for sufficient oxygens (at least 4 for a basic phospho sugar)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, "Too few oxygens for a phospho sugar"

    return True, "Contains a monosaccharide backbone with a phosphate group esterified to a hydroxyl group"