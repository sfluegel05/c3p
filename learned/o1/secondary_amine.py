"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a nitrogen atom bonded to two carbon atoms and zero or one hydrogen,
    excluding cases where nitrogen is part of amides, nitro groups, nitriles, or other
    functional groups where nitrogen is double-bonded to carbon or bonded to oxygen in nitro groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns to exclude unwanted functional groups
    amide_pattern = Chem.MolFromSmarts("N-C(=O)")
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitrile_pattern = Chem.MolFromSmarts("C#N")
    imine_pattern = Chem.MolFromSmarts("N=C")
    azide_pattern = Chem.MolFromSmarts("N=[N+]=[N-]")

    # Check for excluded functional groups
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Contains amide group"
    if mol.HasSubstructMatch(nitro_pattern):
        return False, "Contains nitro group"
    if mol.HasSubstructMatch(nitrile_pattern):
        return False, "Contains nitrile group"
    if mol.HasSubstructMatch(imine_pattern):
        return False, "Contains imine group"
    if mol.HasSubstructMatch(azide_pattern):
        return False, "Contains azide group"

    # Iterate over nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:  # Nitrogen atoms only
            continue

        if atom.GetFormalCharge() != 0:
            continue

        # Check if nitrogen is connected to exactly two carbons
        neighbors = atom.GetNeighbors()
        carbon_count = 0
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 6:
                carbon_count += 1

        if carbon_count != 2:
            continue

        # Check total valence (should be 3 for secondary amine)
        if atom.GetTotalValence() != 3:
            continue

        # Secondary amines can have 0 or 1 hydrogen (due to substitutions like nitrosamines)
        hydrogen_count = atom.GetTotalNumHs()

        if hydrogen_count > 1:
            continue

        # Exclude nitrogens bonded to oxygen (e.g., N-oxides)
        has_oxygen_neighbor = any(nbr.GetAtomicNum() == 8 for nbr in neighbors)
        if has_oxygen_neighbor and not atom.IsInRing():
            continue  # Allow oxygen if in ring (e.g., oximes)

        # If all checks passed, it's a secondary amine
        return True, "Contains secondary amine group"

    return False, "Does not contain secondary amine group"