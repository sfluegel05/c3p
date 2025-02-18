"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    Glucosinolates have a central C atom which is bonded via an S atom to a glycone group
    and via an N atom to a sulfonated oxime group, and which also carries a variable side-group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the glycone pattern (thioglucoside structure)
    thioglucoside_pattern = Chem.MolFromSmarts("S[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(thioglucoside_pattern):
        return False, "No thioglucoside group found"

    # Look for the oxime group connected to sulfonate (R-N=O-S(O)(=O)=O)
    sulfonated_oxime_pattern = Chem.MolFromSmarts("[CX3](=N[OX1][SX4](=[OX1])(=[OX1])[O-])[S]")
    if not mol.HasSubstructMatch(sulfonated_oxime_pattern):
        return False, "No sulfonated oxime group connected"

    # Confirm the central carbon atom with side group
    central_c_pattern = Chem.MolFromSmarts("[CX3](S)~[NX2]")
    if not mol.HasSubstructMatch(central_c_pattern):
        return False, "Central carbon with required connectivity not found"

    # Check stereochemistry for appropriate configuration if possible
    # Implemented stereochemistry checks to match examples could add complexity

    return True, "Molecule contains the defining features of glucosinolate"