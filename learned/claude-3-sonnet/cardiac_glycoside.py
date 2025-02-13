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

    # Check for steroid core
    # More specific steroid pattern with correct connectivity
    steroid_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]([C]1)[C]1[C][C][C]3[C]([C]1[C][C]2)[C][C][C]2[C][C][C][C][C]23")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core found"

    # Check for lactone/butenolide ring (furan-2-one)
    lactone_pattern = Chem.MolFromSmarts("C1=CC(=O)OC1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Check for sugar moiety patterns
    # Look for pyranose ring with multiple OH groups
    sugar_pattern1 = Chem.MolFromSmarts("O[C@H]1[CH2]O[CH]([CH][CH]([CH]1O)O)") # pyranose
    sugar_pattern2 = Chem.MolFromSmarts("O[C@H]1O[CH]([CH]([CH]([CH]1O)O)O)") # alternative sugar pattern
    
    has_sugar = mol.HasSubstructMatch(sugar_pattern1) or mol.HasSubstructMatch(sugar_pattern2)
    if not has_sugar:
        return False, "No sugar moiety found"

    # Check for glycosidic linkage
    glycosidic_pattern = Chem.MolFromSmarts("[C][OH0][C]([CH2,CH])[O][C]") # C-O-C where middle O is not OH
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"

    # Additional checks for typical features
    # Count oxygen atoms (cardiac glycosides typically have many OH groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Too few oxygen atoms for cardiac glycoside"

    # Success - all required features found
    return True, "Contains steroid core with lactone ring, sugar residue(s), and glycosidic linkage"