"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    Glucosinolates are characterized by a central C atom bonded to a glycone group via S,
    a sulfonated oxime group via N, and a side-group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for central carbon bonded to a thioglucoside (sugar) group, which may include stereochemistry
    glycone_pattern = Chem.MolFromSmarts("S[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(glycone_pattern):
        return False, "No thioglucoside group found"
    
    # Look for central carbon bonded to a sulfonated oxime group, specifically N=O-S(O)(=O)=O
    sulfonated_oxime_pattern = Chem.MolFromSmarts("[CX3](=N\\O[S](=O)(=O)[O-])")
    if not mol.HasSubstructMatch(sulfonated_oxime_pattern):
        return False, "No sulfonated oxime group connected"
    
    # Confirm central carbon atom connectivity pattern
    # Bonded via sulfur to glycone and via nitrogen to sulfonated oxime, additionally carrying a side chain
    central_c_pattern = Chem.MolFromSmarts("[CX3](S)(N=O)~[*]")
    if not mol.HasSubstructMatch(central_c_pattern):
        return False, "Central carbon with required connectivity not found"
    
    return True, "Molecule contains the defining features of glucosinolate"