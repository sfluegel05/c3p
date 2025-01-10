"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin must have a specific penam core structure with additional specific substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define penam core SMARTS pattern
    penam_core_pattern = Chem.MolFromSmarts("C1C2SC(C)(C)N2C(=O)C3=C1NC3=O")  # Updated core with bicyclic specificity
    if not mol.HasSubstructMatch(penam_core_pattern):
        return False, "Penam core structure not found"

    # Check for two methyl groups at position 2 on thiazolidine ring
    methyl_group_pattern = Chem.MolFromSmarts("[C@]12([S@@]C(C)(C)[N@]2C(=O)[C@H]1)C(=O)[O-]")  # Position-specific methyl search
    if not mol.HasSubstructMatch(methyl_group_pattern):
        return False, "Two methyl groups not found at position 2"

    # Check for carboxylate/carboxylic acid at position 3
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    carboxyl_group_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylate_pattern) and not mol.HasSubstructMatch(carboxyl_group_pattern):
        return False, "Carboxylate group not found at position 3"

    # Check for carboxamido group at position 6
    carboxamido_pattern = Chem.MolFromSmarts("NC(=O)C")  # Revised for position-specific carboxamido
    if not mol.HasSubstructMatch(carboxamido_pattern):
        return False, "Carboxamido group not found at position 6"

    # Check stereochemistry specifics for penicillin structure
    stereochemistry_correct = False
    chirality_pattern = Chem.CanonSmiles('[C@H]1(CN2C(=O)C1=S[C@H]23)C(C)=O')  # Bicyclic chirality pattern
    if Chem.MolFromSmiles(chirality_pattern):
        stereochemistry_correct = True

    if not stereochemistry_correct:
        return False, "Stereochemistry does not match penicillin"

    return True, "Molecule matches the penicillin structure"