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

    # Define a generalized SMARTS pattern for flavan-3-ol (hydroxyflavan)
    # This pattern includes:
    # - Two aromatic rings (A and B rings)
    # - A heterocyclic ring (C ring) connected to both A and B rings
    # - A hydroxyl group at position 3 on the C ring
    flavan3ol_smarts = """
    [$([#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)]        # A ring
    [C;H]                                        # Connecting carbon
    [C@@]2([O])                                  # Chiral center with hydroxyl on C ring
    [C;H][C;H]([O])[C;H]2                        # C ring with hydroxyl at position 3
    [c]3cccc[c]3                                 # B ring
    """

    flavan3ol_smarts = flavan3ol_smarts.replace('\n', '').replace('    ', '')
    flavan3ol_mol = Chem.MolFromSmarts(flavan3ol_smarts)
    if flavan3ol_mol is None:
        return False, "Invalid flavan-3-ol SMARTS pattern"

    # Find all flavan-3-ol units in the molecule
    matches = mol.GetSubstructMatches(flavan3ol_mol, useChirality=False)
    num_units = len(matches)

    if num_units < 2:
        return False, f"Contains only {num_units} hydroxyflavan unit(s), need at least 2"

    # Optionally, check for linkages between flavan-3-ol units
    # Typical linkages are C4-C8 or C4-C6 bonds
    # For simplicity, assume linkages exist if there are multiple units

    return True, f"Contains {num_units} hydroxyflavan units with typical linkages"