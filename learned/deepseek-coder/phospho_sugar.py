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

    # Look for phosphate group pattern (P with 3 or 4 oxygens, one of which is double-bonded)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group found"

    # Look for sugar backbone (3-8 carbons with multiple hydroxyl groups)
    # More flexible pattern that matches both linear and cyclic sugars
    sugar_pattern = Chem.MolFromSmarts("[C;H0,H1,H2][C;H0,H1,H2]([OH])[C;H0,H1,H2]([OH])")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        # Check for ring structures with at least two hydroxyl groups
        ring_sugar_pattern = Chem.MolFromSmarts("[C;H0,H1,H2]1[C;H0,H1,H2][C;H0,H1,H2]([OH])[C;H0,H1,H2]([OH])[C;H0,H1,H2]1")
        ring_sugar_matches = mol.GetSubstructMatches(ring_sugar_pattern)
        if len(ring_sugar_matches) == 0:
            return False, "No sugar backbone found"

    # Check if the phosphate is attached to the sugar via an ester bond
    # More comprehensive ester bond pattern
    ester_pattern = Chem.MolFromSmarts("[C;H0,H1,H2][OX2][PX4](=[OX1])([OX2])[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "Phosphate not esterified to sugar"

    # Count carbons and oxygens to ensure it's a monosaccharide
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Relax carbon count restriction to 3-8 to accommodate heptoses
    if c_count < 3 or c_count > 8:
        return False, "Not a monosaccharide (incorrect number of carbons)"
    if o_count < 4:
        return False, "Too few oxygens for a phospho sugar"

    return True, "Contains a monosaccharide backbone with a phosphate group esterified to a hydroxyl group"