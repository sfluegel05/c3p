"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Define SMARTS pattern for flavan-3-ol (catechin/epicatechin) unit
    # This pattern matches the core structure of hydroxyflavans
    flavan3ol_smarts = """
    [#6]1=CC(=CC=C1)[C@@H]2O[C@H](C[C@H](O)[C@@H]2O)c3cc(O)ccc3
    """

    flavan3ol_mol = Chem.MolFromSmarts(flavan3ol_smarts)
    if flavan3ol_mol is None:
        return False, "Invalid flavan-3-ol SMARTS pattern"

    # Find all flavan-3-ol units in the molecule
    matches = mol.GetSubstructMatches(flavan3ol_mol, useChirality=True)
    num_units = len(matches)

    if num_units >= 2:
        return True, f"Contains {num_units} hydroxyflavan (flavan-3-ol) units"
    else:
        return False, f"Contains only {num_units} hydroxyflavan unit(s), need at least 2"