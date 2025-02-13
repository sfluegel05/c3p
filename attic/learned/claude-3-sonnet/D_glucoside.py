"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: D-glucoside
A glucoside where the glycoside group is derived from D-glucose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for beta-D-glucose core with correct stereochemistry
    # [C@@H]-O-[C@H]-[C@@H]-[C@H]-[C@H]-CH2OH with correct OH positions
    glucose_pattern = Chem.MolFromSmarts('[C@@H]1([OR0])[C@H]([OR0])[C@@H]([OR0])[C@H]([OR0])[C@H](CO[R0])O1')
    
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No D-glucose core found with correct stereochemistry"
    
    # Check for glycosidic bond (O-linkage)
    glycosidic_pattern = Chem.MolFromSmarts('[CR1]-O-[CR1]')
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic bond found"
    
    # Count OH groups on glucose ring (should have multiple)
    oh_pattern = Chem.MolFromSmarts('[OH1]')
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    
    if oh_matches < 3:  # Glucose should have multiple OH groups
        return False, "Too few hydroxyl groups for a glucoside"
        
    # Additional check for cyclic structure
    ring_pattern = Chem.MolFromSmarts('[C]1[C][C][C][C][C]1')
    if not mol.HasSubstructMatch(ring_pattern):
        return False, "No six-membered ring found"
    
    return True, "Contains D-glucose core with glycosidic linkage and correct stereochemistry"