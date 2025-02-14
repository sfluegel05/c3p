"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid with an ester linkage formed by condensation of a carboxylic acid
    with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General steroid (four ring pattern) with connections
    steroid_pattern = Chem.MolFromSmarts("C1CC[C@H]2[C@@H]3CC[C@H]4[C@@]3(CC[C@@H]2[C@@H]1C4)C")
    
    # Ester pattern - specifically connected to an oxygen (3-hydroxy position potential)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)")
    
    # Check for steroid skeleton
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid skeleton found"

    # Check for ester linkage at specific position (3-hydroxy)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "Ester linkage not found at the 3-hydroxy group"

    # If both steroid skeleton and correct ester linkage were found
    return True, "Contains a steroid skeleton with the ester linkage on the 3-hydroxy group"

# Test the function with a sterol ester example
smiles_example = "CCCCCCCCCCCCCCCCC(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)[C@H](C)CCCC(C)C"
print(is_sterol_ester(smiles_example))