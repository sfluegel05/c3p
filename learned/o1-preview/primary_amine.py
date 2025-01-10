"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine has a nitrogen atom bonded to two hydrogens and one carbon (NH2 group attached to a hydrocarbyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to the molecule
    mol = Chem.AddHs(mol)

    # Define primary amine pattern: nitrogen with two hydrogens and single bond to carbon
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2][CX4]")

    # Define amide pattern: nitrogen connected to a carbonyl carbon
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")

    # Define nitro group pattern
    nitro_pattern = Chem.MolFromSmarts("[NX3](=O)[O-]")

    # Find amide nitrogens
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    amide_nitrogens = set()
    for match in amide_matches:
        nitrogen_idx = match[0]
        amide_nitrogens.add(nitrogen_idx)

    # Find nitro nitrogens
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    nitro_nitrogens = set()
    for match in nitro_matches:
        nitrogen_idx = match[0]
        nitro_nitrogens.add(nitrogen_idx)

    # Search for primary amine groups
    primary_amines = mol.GetSubstructMatches(primary_amine_pattern)
    if not primary_amines:
        return False, "No primary amine group found"

    # Check each primary amine to ensure it's not part of an amide or nitro group
    for match in primary_amines:
        nitrogen_idx = match[0]
        if nitrogen_idx in amide_nitrogens:
            continue  # Skip amide nitrogen
        if nitrogen_idx in nitro_nitrogens:
            continue  # Skip nitro nitrogen

        # Ensure nitrogen is not in a ring (exclude heterocyclic amines)
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        if nitrogen_atom.IsInRing():
            continue

        # If all checks passed, it's a primary amine
        return True, "Contains primary amine group (NH2 attached to hydrocarbyl group)"

    return False, "No valid primary amine group found"

__metadata__ = {
    'chemical_class': {
        'name': 'primary amine',
        'definition': 'A compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.'
    },
    'config': {
        'additional_info': 'Classifies primary amines based on NH2 group attached to hydrocarbyl group.'
    }
}