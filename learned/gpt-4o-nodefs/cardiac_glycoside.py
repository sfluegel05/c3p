"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid core recognition (cyclopenta[a]phenanthrene skeleton)
    steroid_patterns = [
        Chem.MolFromSmarts("C1CCC2C3CCC4C(C3CC2C1)CCC4"),
        Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCC4")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return False, "No steroid core found"

    # Glycoside detection (O-glycosidic linkage)
    glycoside_pattern = Chem.MolFromSmarts("C-O[C@@H]1[C@H]([C@H]([C@@H]([C@@H](O1)CO)O)O)O")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosidic linkage found"

    # Lactone ring detection
    lactone_patterns = [
        Chem.MolFromSmarts("O=C1COC=C1"),
        Chem.MolFromSmarts("O=C1COCC1")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in lactone_patterns):
        return False, "No lactone ring found"
    
    # Hydroxyl groups verification
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetTotalDegree() == 1)
    if hydroxyl_count < 3:
        return False, f"Insufficient hydroxyl groups in the molecule: {hydroxyl_count} found"
    
    return True, "The structure contains a steroid core with glycosidic linkage and a lactone ring characteristic of cardiac glycosides"