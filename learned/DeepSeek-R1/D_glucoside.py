"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside has a beta-D-glucopyranose group linked via a glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # SMARTS pattern for beta-D-glucopyranose with glycosidic bond
    # Matches the glucopyranose ring with O-linked group at C1
    glucoside_pattern = Chem.MolFromSmarts("[C@H]1([C@@H](O)[C@H](O)[C@@H](O)[C@@H](O1)CO)O[!H]")
    if mol.HasSubstructMatch(glucoside_pattern):
        return True, "Contains beta-D-glucopyranose with glycosidic bond"
    
    # Alternative pattern considering different stereochemistry possibilities
    # Check for any glucopyranose with O-linkage at C1 and correct D-configuration
    # This pattern may need adjustment based on specific examples
    glucoside_pattern2 = Chem.MolFromSmarts("[OX2][C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O[!H]")
    if mol.HasSubstructMatch(glucoside_pattern2):
        return True, "Contains beta-D-glucopyranose with glycosidic bond"
    
    return False, "No beta-D-glucopyranose with glycosidic bond found"