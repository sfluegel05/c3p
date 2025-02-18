"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmarts

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine has a benzene ring with two adjacent hydroxyl groups (catechol)
    and a 2-aminoethyl group (-CH2-CH2-NR2) attached to the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for benzene ring with two adjacent hydroxyl groups (catechol)
    catechol_pattern = MolFromSmarts("[OH]c1c([OH])cccc1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol group (benzene-1,2-diol) found"

    # Check for 2-aminoethyl group (-CH2-CH2-NR2) attached to aromatic carbon
    # The pattern allows substituents on the ethyl chain and nitrogen
    ethylamine_pattern = MolFromSmarts("[c](-[CH2]-[CH2]-[NX3])")
    if not mol.HasSubstructMatch(ethylamine_pattern):
        return False, "No 2-aminoethyl group attached to aromatic ring"

    return True, "Contains catechol group with 2-aminoethyl substituent"