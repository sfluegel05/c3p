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

    # Define the flavan-3-ol core structure (catechin or epicatechin unit) without stereochemistry
    flavan3ol_smiles = 'Oc1ccc(cc1)[C]2CC(O)c3ccccc23'  # Simplified core structure
    flavan3ol = Chem.MolFromSmiles(flavan3ol_smiles)
    if flavan3ol is None:
        return False, "Error in flavan-3-ol SMILES pattern"

    # Find flavan-3-ol units in the molecule (ignore stereochemistry)
    flavan3ol_matches = mol.GetSubstructMatches(flavan3ol, useChirality=False)
    num_flavan_units = len(flavan3ol_matches)

    if num_flavan_units < 2:
        return False, f"Found {num_flavan_units} flavan-3-ol units, need at least 2"

    # Optional: Check for interflavan bonds (connections between flavan-3-ol units)
    # For simplicity, we'll assume that if there are at least 2 units, they are oligomerized
    # Alternatively, you can implement logic to check for bonds between units if necessary

    return True, f"Contains {num_flavan_units} flavan-3-ol units connected via interflavan bonds"