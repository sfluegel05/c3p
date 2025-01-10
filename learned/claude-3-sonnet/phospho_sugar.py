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

    # Look for phosphate group pattern
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[O])([OX2H,OX1-])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for sugar patterns - both cyclic and open chain forms
    # Furanose pattern (5-membered ring with oxygens)
    furanose_pattern = Chem.MolFromSmarts("[CH2X4,CH1X4,CH0X4]-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[OX2]")
    # Pyranose pattern (6-membered ring with oxygens)
    pyranose_pattern = Chem.MolFromSmarts("[CH2X4,CH1X4,CH0X4]-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[OX2]")
    # Open chain sugar pattern (multiple hydroxyls)
    open_sugar_pattern = Chem.MolFromSmarts("[CH2X4,CH1X4](-[OX2])-[CH1X4](-[OX2])-[CH1X4](-[OX2])-[CH1X4](-[OX2])")

    is_sugar = False
    if mol.HasSubstructMatch(furanose_pattern):
        is_sugar = True
    elif mol.HasSubstructMatch(pyranose_pattern):
        is_sugar = True
    elif mol.HasSubstructMatch(open_sugar_pattern):
        is_sugar = True

    if not is_sugar:
        return False, "No sugar moiety found"

    # Verify phosphate is connected to sugar via ester linkage
    phospho_ester_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[P](=[O])([OX2H,OX1-])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(phospho_ester_pattern):
        return False, "Phosphate not connected to sugar via ester bond"

    # Count hydroxyl groups (excluding phosphate hydroxyls)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4]-[OX2H]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches < 2:
        return False, "Insufficient hydroxyl groups for a sugar"

    # Additional check for nucleotides - if has nucleobase, must have sugar-phosphate
    nucleobase_pattern = Chem.MolFromSmarts("[#7r5,#7r6]1[#6r5,#6r6][#7r5,#7r6][#6r5,#6r6][#6r5,#6r6]1")
    if mol.HasSubstructMatch(nucleobase_pattern):
        sugar_phosphate_pattern = Chem.MolFromSmarts("[CH2X4]-[OX2]-[P](=[O])([OX2H,OX1-])[OX2H,OX1-]")
        if not mol.HasSubstructMatch(sugar_phosphate_pattern):
            return False, "Nucleotide without proper sugar-phosphate linkage"

    return True, "Contains monosaccharide with phosphate ester"