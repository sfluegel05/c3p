"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    Glucosinolates have a thioglucoside linkage, a central C linked via S
    and N to a sulfonated oxime group, and carry a side-group.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a glucosinolate, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized thioglucoside linkage pattern 
    # (more flexible with stereocenters and linkage possibilities)
    thioglucoside_pattern = Chem.MolFromSmarts("S[C@H1]1OC(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(thioglucoside_pattern):
        return False, "No thioglucoside linkage found"

    # Match the sulfonated oxime group
    sulfonated_oxime_pattern = Chem.MolFromSmarts("C=N/OS(=O)(=O)[O-]")
    if not mol.HasSubstructMatch(sulfonated_oxime_pattern):
        return False, "Sulfonated oxime linkage not found"

    # Ensure the central C is linked via S and N, and extends to a side-group
    central_pattern = Chem.MolFromSmarts("S[C]=N/OS(=O)(=O)[O-]")
    if not mol.HasSubstructMatch(central_pattern):
        return False, "Central carbon linkage with S and N not correctly represented"

    return True, "Contains all structural features of a glucosinolate"