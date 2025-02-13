"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides have a steroid core with a lactone ring and sugar residues.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid core - more flexible pattern
    # Uses cyclopentanoperhydrophenanthrene core with flexible connectivity
    steroid_pattern = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~3~4~2")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core found"

    # Check for butenolide (unsaturated lactone) ring - characteristic of cardiac glycosides
    lactone_pattern = Chem.MolFromSmarts("O=C1OCC=C1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No butenolide ring found"

    # Check for sugar moieties - multiple patterns to catch different sugar types
    sugar_patterns = [
        Chem.MolFromSmarts("O[C@H]1[CH2]O[CH]([CH][CH]([CH]1O)O)"), # pyranose
        Chem.MolFromSmarts("O[C@H]1O[CH]([CH]([CH]([CH]1O)O)O)"),   # alternative pyranose
        Chem.MolFromSmarts("OC1C(O)C(O)C(O)CO1"),                    # basic sugar ring
        Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1")                 # fully hydroxylated sugar
    ]
    
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_sugar:
        return False, "No sugar moiety found"

    # Check for glycosidic linkage at position 3 of steroid
    glycosidic_patterns = [
        Chem.MolFromSmarts("[C]1[C][C]([OH0][C])[C]"), # basic glycosidic bond
        Chem.MolFromSmarts("[C]1[C][C]([OH0][CH]2O[CH])[C]") # more specific pattern
    ]
    
    has_glycosidic = any(mol.HasSubstructMatch(pattern) for pattern in glycosidic_patterns)
    if not has_glycosidic:
        return False, "No glycosidic linkage found"

    # Additional structural features common in cardiac glycosides
    
    # Check for hydroxyl groups - typically at positions 3, 14, and others
    oh_pattern = Chem.MolFromSmarts("[C]([C][C])([C][C])[OH]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 2:
        return False, "Insufficient hydroxyl groups for cardiac glycoside"

    # Count oxygen atoms (cardiac glycosides typically have many oxygen atoms)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 8:  # increased from 6 to be more specific
        return False, "Too few oxygen atoms for cardiac glycoside"

    # Check molecular weight - cardiac glycosides are typically large molecules
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 400:  # typical cardiac glycosides are >500 Da
        return False, "Molecular weight too low for cardiac glycoside"

    # Success - all required features found
    return True, "Contains steroid core with butenolide ring, sugar residue(s), and characteristic glycosidic linkage"