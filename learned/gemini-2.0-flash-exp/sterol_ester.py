"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid with an ester bond at the 3-hydroxy group.

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

    # 1. Check for steroid core with ester at position 3
    # This pattern matches the core fused rings, with possible substitutions,
    # and an ester bond connected to the carbon at position 3.
    # we are using the "a" symbol instead of "C" to allow for double bonds.
    sterol_ester_pattern = Chem.MolFromSmarts("[a]1[a][a][CX4]([OX2][CX3](=[OX1])[#6])[a]2[a][a][a]3[a][a][a][a]4[a][a]3[a]214")
    if not mol.HasSubstructMatch(sterol_ester_pattern):
        return False, "No steroid core with ester at position 3 found"


    return True, "Sterol ester detected: Contains a steroid core with an ester bond at position 3"