"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    Glucosinolates are defined by a central C atom bonded to a glycone group via S,
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

    # Glycone thioglucoside group pattern (allow stereochemical flexibility)
    glycone_pattern = Chem.MolFromSmarts("S[C@@H]1O[C@@H]([C@@H](O)[C@@H]1O)CO")
    if not mol.HasSubstructMatch(glycone_pattern):
        return False, "No thioglucoside group found"

    # N-sulfooxy group pattern :
    sulfonated_oxime_pattern = Chem.MolFromSmarts("[CX3](=N\\OS(=O)(=O)[O-])")
    oxime_matches = mol.GetSubstructMatches(sulfonated_oxime_pattern)
    if len(oxime_matches) == 0:
        return False, "No sulfonated oxime group connected"

    # Central carbon connectivity, ensuring linkage to both S and N features
    central_c_pattern = Chem.MolFromSmarts("C(S)(=N\\O[S](=O)(=O)[O-])")
    if not mol.HasSubstructMatch(central_c_pattern):
        return False, "Central carbon with required connectivity not found"

    return True, "Molecule contains the defining features of glucosinolate"