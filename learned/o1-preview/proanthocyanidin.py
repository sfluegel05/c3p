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

    # Define SMARTS pattern for flavan-3-ol (without stereochemistry)
    # This pattern matches the core structure of hydroxyflavans
    flavan3ol_smarts = 'c1cc(O)ccc1[C]2OC[C](O)C2'  # Simplified flavan-3-ol core

    flavan3ol_mol = Chem.MolFromSmarts(flavan3ol_smarts)
    if flavan3ol_mol is None:
        return False, "Invalid flavan-3-ol SMARTS pattern"

    # Find all flavan-3-ol units in the molecule
    matches = mol.GetSubstructMatches(flavan3ol_mol, useChirality=False)
    num_units = len(matches)

    if num_units < 2:
        return False, f"Contains only {num_units} hydroxyflavan unit(s), need at least 2"

    # Now check for linkages between flavan-3-ol units
    # Typical linkages are C4-C8 or C4-C6 bonds
    # Define possible linkage patterns
    linkage_patterns = [
        '[C;R1]1([C;R1])[C;R1][C;R1][C;R1][C;R1][C;R1]1',  # Ring structure
        '[C;R1][C;R1][C;R1][C;R1][C;R1][C;R1]'             # Chain structure
    ]

    has_linkage = False
    for patt in linkage_patterns:
        linkage_mol = Chem.MolFromSmarts(patt)
        if linkage_mol is not None and mol.HasSubstructMatch(linkage_mol):
            has_linkage = True
            break

    if not has_linkage:
        return False, "No typical linkage between flavan-3-ol units found"

    return True, f"Contains {num_units} hydroxyflavan units with typical linkages"