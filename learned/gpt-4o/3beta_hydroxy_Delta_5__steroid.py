"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid is a 3beta-hydroxy-steroid that contains a double bond
    between positions 5 and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 3beta-hydroxy group SMARTS pattern, aiming for general steroids:
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H]([C@H2][C@H](C)[C@H]=O)(O)[C@@H](C)C")

    # Steroid nucleus criterion: four fused rings (ABCD rings, common in steroids)
    steroid_pattern = Chem.MolFromSmarts(
        "[#6]1([#6][#6][#6]2[#6][#6][#8][#6]3[#6][#6][#8][#6]4[#6][#6][#6][#6][#6][#6][#6]4[C@H]3[C@@H]1)"
    )

    # Delta(5) specific double bond SMARTS pattern:
    delta5_pattern = Chem.MolFromSmarts("[C;!R]=[C;!R&$(C1C(C=CC(C)=C1)=O)]")

    # Check for 3beta-hydroxy group
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxy group match found."

    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone detected."

    # Check for Delta(5) double bond
    if not mol.HasSubstructMatch(delta5_pattern):
        return False, "No specific Delta(5) double bond match found."

    return True, "Molecule is classified as 3beta-hydroxy-Delta(5)-steroid."