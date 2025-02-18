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
    and a 2-aminoethyl group (-CH2-CH2-NR2) or its substituted derivatives attached to the ring.

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
    # SMARTS: two adjacent aromatic hydroxyls in a benzene ring
    catechol_pattern = MolFromSmarts("[OH][c;r6]@[c;r6][OH]")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol group (adjacent aromatic hydroxyls) found"

    # Check for 2-aminoethyl group or substituted derivatives
    # SMARTS: aromatic carbon connected to two carbons (any substitution) followed by an amine (not amide)
    ethylamine_pattern = MolFromSmarts("[c;r6]-[CX4]-[CX4]-[NX3;H0,H1,H2;!$(NC=O)]")
    if not mol.HasSubstructMatch(ethylamine_pattern):
        return False, "No 2-aminoethyl-like group attached to aromatic ring"

    return True, "Contains catechol group with 2-aminoethyl-like substituent"