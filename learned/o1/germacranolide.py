"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on a germacrane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for germacrane skeleton (10-membered ring with specific connectivity)
    germacrane_smarts = "[C;R]1[C;R][C;R][C;R][C;R][C;R][C;R][C;R][C;R][C;R]1"  # 10-membered ring
    germacrane_mol = Chem.MolFromSmarts(germacrane_smarts)
    if not mol.HasSubstructMatch(germacrane_mol):
        return False, "No germacrane skeleton (10-membered ring) found"

    # Define SMARTS pattern for lactone ring attached to germacrane skeleton
    lactone_smarts = "[C;R](=O)[O;R][C;R]"  # Lactone group within a ring
    lactone_mol = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_mol):
        return False, "No lactone group (cyclic ester) found within the ring"

    # Ensure lactone is part of the 10-membered ring
    germacrane_matches = mol.GetSubstructMatch(germacrane_mol)
    lactone_matches = mol.GetSubstructMatches(lactone_mol)
    lactone_in_ring = False
    for lactone_match in lactone_matches:
        if set(lactone_match).issubset(germacrane_matches):
            lactone_in_ring = True
            break
    if not lactone_in_ring:
        return False, "Lactone group is not part of the germacrane ring"

    # Optionally, check for sesquiterpene backbone (15 carbons) but allow for substitutions
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Too few carbons ({c_count}) for a germacranolide (sesquiterpene)"

    return True, "Molecule is a germacranolide (germacrane skeleton with lactone ring)"