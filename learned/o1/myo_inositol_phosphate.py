"""
Classifies: CHEBI:25448 myo-inositol phosphate
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

    # Define the SMARTS pattern for the myo-inositol core with variable substituents on oxygen atoms
    # The pattern includes the chiral centers (useChirality=True in substructure search)
    myo_inositol_phosphate_smarts = "[C@H]1(O[*])[C@@H](O[*])[C@H](O[*])[C@@H](O[*])[C@H](O[*])[C@@H]1(O[*])"
    myo_inositol_query = Chem.MolFromSmarts(myo_inositol_phosphate_smarts)
    if myo_inositol_query is None:
        return None, "Cannot create myo-inositol phosphate query"

    # Perform a substructure match considering stereochemistry
    matches = mol.GetSubstructMatches(myo_inositol_query, useChirality=True)
    if matches:
        return True, "Contains myo-inositol phosphate core structure"
    else:
        return False, "Does not contain myo-inositol phosphate core structure"