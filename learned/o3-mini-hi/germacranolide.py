"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: Germacranolide – a sesquiterpene lactone based on the germacrane skeleton
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    Germacranolides are sesquiterpene lactones that have (1) a lactone ring (cyclic ester)
    and (2) a germacrane skeleton, defined by a 10-membered (mostly carbon) ring.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a germacranolide, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, check for a lactone ring.
    # We define a lactone as a cyclic ester, so we use a SMARTS that finds an ester group
    # (C(=O)O) where both the carbonyl carbon and the oxygen are in a ring.
    lactone_smarts = "[#6;R](=O)[O;R]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring (cyclic ester) found"

    # Next, look for a 10-membered ring characteristic of the germacrane skeleton.
    # We inspect each ring in the molecule to see if there is one with 10 atoms and,
    # further, if at least 8 of those atoms are carbon (allowing for one or two heteroatoms).
    ring_info = mol.GetRingInfo()
    germacrane_found = False
    for ring in ring_info.AtomRings():
        if len(ring) == 10:
            carbon_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if carbon_count >= 8:
                germacrane_found = True
                break
    if not germacrane_found:
        return False, "No 10-membered (germacrane-type) ring found"

    # If both features (lactone ring and 10-membered ring) are present, we classify as germacranolide.
    return True, "Found both a lactone ring and a 10-membered germacrane skeleton typical for germacranolide"