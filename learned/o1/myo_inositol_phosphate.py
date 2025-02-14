"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component has myo-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for myo-inositol core with correct stereochemistry
    # Allow any substituent attached to the oxygen atoms
    myo_inositol_smarts = '[C@@H]1([O][*])[C@H]([O][*])[C@@H]([O][*])[C@H]([O][*])[C@@H]([O][*])[C@@H]1[O][*]'
    pattern = Chem.MolFromSmarts(myo_inositol_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Perform substructure search with chirality consideration
    matches = mol.GetSubstructMatches(pattern, useChirality=True)

    if matches:
        return True, "Contains myo-inositol phosphate core structure"
    else:
        return False, "No myo-inositol phosphate core structure found"