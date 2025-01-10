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

    # Generalized thioglucoside linkage pattern: sugar structure linked through sulfur
    thioglucoside_pattern = Chem.MolFromSmarts("S[C@H]1OC[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(thioglucoside_pattern):
        return False, "No thioglucoside linkage (linked through sulfur) found"

    # Match the sulfonated oxime group: C=N-OS(=O)(=O)[O-], with anti configuration
    sulfonated_oxime_pattern = Chem.MolFromSmarts("C=N/OS(=O)(=O)[O-]")
    if not mol.HasSubstructMatch(sulfonated_oxime_pattern):
        return False, "No appropriate sulfonated oxime linkage found"

    # Match central linking structure with side-chain: include central C-S and C=N 
    central_pattern = Chem.MolFromSmarts("S[C](=N/OS(=O)(=O)[O-])[CX4,CX3]")
    if not mol.HasSubstructMatch(central_pattern):
        return False, "Central carbon structure not recognized with correct linkages and side-chain"
    
    return True, "Contains all structural features of a glucosinolate"