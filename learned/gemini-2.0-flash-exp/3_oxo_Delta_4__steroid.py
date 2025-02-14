"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid has a 3-ketone group and a C=C double bond between position 4 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for steroid core
    # This pattern matches a basic 4-ring system with R-groups at each position
    # Three 6-membered rings and one 5-membered ring
    steroid_core_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CR1][CR1]2[CR1]1[CR1]([CR1])([CR1])[CR1]3[CR1]2[CR1][CR1][CR1]4[CR1]3[CR1][CR1]4")

    # Define SMARTS for the 3-oxo group with delta-4 double bond
    oxo_group_pattern = Chem.MolFromSmarts("[C]1-[C](=[O])-[C]=[C]-[C]1")

    # Check for substructure matches
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"

    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo-Delta(4) group found"

    return True, "Contains steroid core with 3-oxo and Delta(4) double bond."