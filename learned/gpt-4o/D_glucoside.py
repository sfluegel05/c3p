"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is a glucoside in which the glycoside group is derived from D-glucose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define corrected SMARTS for a glucose moiety
    # Using a more generalized depiction of the glucose ring backbone
    d_glucose_smarts = "OC[C@H]1O[C@@H](CO)[C@H](O)[C@H](O)[C@@H]1O"
    d_glucose_mol = Chem.MolFromSmarts(d_glucose_smarts)
    if d_glucose_mol is None:
        return None, "Error in setting up D-glucose SMARTS pattern."
    
    # Check if D-glucose moiety is present
    if not mol.HasSubstructMatch(d_glucose_mol):
        return False, "D-glucose moiety not found"
    
    # Simple check for glycosidic linkage using an anomeric carbon connected to an ether
    # Anomeric carbon will typically be attached to another oxygen and further ether bonds 
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 4:
            oxygen_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 8]
            if len(oxygen_neighbors) == 2:  # Connected to two oxygens for glycosidic bond
                # Assuming one is within the ring and another outside the ring (glycosidic bond)
                return True, "Contains D-glucose moiety with glycosidic linkage"

    return False, "Glycosidic linkage not found to D-glucose moiety"