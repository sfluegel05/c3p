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
    
    # Define SMARTS for D-glucose moiety having either alpha or beta configuration
    d_glucose_smarts = '[C@H]1(O)[C@@H]([C@H](O)[C@@H](CO)O)O[C@@H]1O |#1:3@1:0,@3:6&@4|'
    d_glucose_mol = Chem.MolFromSmarts(d_glucose_smarts)
    
    # Check if D-glucose moiety is present
    if not mol.HasSubstructMatch(d_glucose_mol):
        return False, "D-glucose moiety not found"
    
    # Look for any glycosidic linkage by checking for ether bond from anomeric carbon
    anomeric_c = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 3:
            connected_oxygen = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    connected_oxygen = True
            if connected_oxygen:
                anomeric_c = atom
                break

    if anomeric_c is None:
        return False, "Anomeric carbon not found; probably not a glycoside"

    # Check for ether linkage from anomeric carbon
    for neighbor in anomeric_c.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:
            if len([n for n in neighbor.GetNeighbors() if n.GetIdx() != anomeric_c.GetIdx()]) > 0:
                return True, "Contains D-glucose moiety with glycosidic linkage"

    return False, "Glycosidic linkage not found to D-glucose moiety"