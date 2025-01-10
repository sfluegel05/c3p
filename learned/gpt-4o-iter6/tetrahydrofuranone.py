"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is an oxolane with an oxo- substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # The SMARTS pattern for tetrahydrofuranone
    # The pattern specifies a 5-membered ring with one oxygen and one carbonyl group
    triangle_oxygen_carbonyl = Chem.MolFromSmarts("C1(=O)OC(C)C1")
    generic_tetrahydrofuranone_pattern = Chem.MolFromSmarts("C1(=O)OCC[C@1]")
    
    # Match the structure
    if mol.HasSubstructMatch(triangle_oxygen_carbonyl) or mol.HasSubstructMatch(generic_tetrahydrofuranone_pattern):
        return True, "Tetrahydrofuranone structure confirmed"

    return False, "No valid tetrahydrofuranone structure found"