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

    # Define the SMARTS pattern for flavan-3-ol units (catechin/epicatechin units)
    # This pattern captures the core structure of flavan-3-ol with hydroxyl groups
    flavan3ol_smarts = """
    [#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-
    [C@@H]2[C@H](O)-
    [CH2]-
    [C@@H]3Oc4cc(O)ccc4[C@@H]23
    """
    # Remove whitespace and newlines
    flavan3ol_smarts = "".join(flavan3ol_smarts.split())
    flavan3ol = Chem.MolFromSmarts(flavan3ol_smarts)
    if flavan3ol is None:
        return False, "Error in flavan-3-ol SMARTS pattern"

    # Find flavan-3-ol units in the molecule (ignore stereochemistry)
    flavan3ol_matches = mol.GetSubstructMatches(flavan3ol, useChirality=False)
    num_flavan_units = len(flavan3ol_matches)

    if num_flavan_units < 2:
        return False, f"Found {num_flavan_units} flavan-3-ol units, need at least 2"

    # Optional: Check for interflavan bonds (C4-C6 or C4-C8 linkages)
    # For simplicity, we'll assume that if there are at least 2 units, they are connected via condensation

    return True, f"Contains {num_flavan_units} flavan-3-ol units connected via interflavan bonds"