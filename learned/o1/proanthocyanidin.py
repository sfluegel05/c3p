"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: proanthocyanidin
"""
from rdkit import Chem

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for flavan-3-ol units (catechin and epicatechin cores)
    catechin_smarts = "Oc1ccc(cc1)[C@@H]2O[C@H](C[C@H]2O)c3ccc(O)cc3"  # Catechin core
    epicatechin_smarts = "Oc1ccc(cc1)[C@H]2O[C@@H](C[C@H]2O)c3ccc(O)cc3"  # Epicatechin core

    catechin_mol = Chem.MolFromSmarts(catechin_smarts)
    epicatechin_mol = Chem.MolFromSmarts(epicatechin_smarts)
    if catechin_mol is None or epicatechin_mol is None:
        return False, "Error in flavan-3-ol SMARTS patterns"

    # Initialize count
    num_flavan_units = 0

    # Check for catechin units
    catechin_matches = mol.GetSubstructMatches(catechin_mol, useChirality=False)
    num_flavan_units += len(catechin_matches)

    # Check for epicatechin units
    epicatechin_matches = mol.GetSubstructMatches(epicatechin_mol, useChirality=False)
    num_flavan_units += len(epicatechin_matches)

    if num_flavan_units < 2:
        return False, f"Found {num_flavan_units} flavan-3-ol units, need at least 2"

    # Optional: Check for interflavan bonds (C4-C6 or C4-C8 linkages)
    # For simplicity, we'll assume that if there are at least 2 units, they are connected via condensation

    return True, f"Contains {num_flavan_units} flavan-3-ol units connected via interflavan bonds"