"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string provided."

    # Define the SMARTS pattern for the butenolide core (2-furanone)
    # This matches 5-membered rings with one oxygen, a carbonyl group, and adjacent double-bonded carbon
    butenolide_pattern = Chem.MolFromSmarts("[OX2]1[CX3](=[OX1])[CX3]=[CX3]1")

    # Check if the molecule has the butenolide core structure
    if mol.HasSubstructMatch(butenolide_pattern):
        return True, "Molecule contains the 2-furanone core structure."
    else:
        return False, "Molecule does not contain the 2-furanone core structure."