"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: CHEBI: monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine contains one amino group connected to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find aromatic rings
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}
    if not aromatic_atoms:
        return False, "No aromatic ring found"

    # SMARTS pattern for amino group connected via two-carbon chain to aromatic ring
    # Pattern: aromatic carbon -> any carbon -> nitrogen with at least one H (primary, secondary amine)
    # Excludes tertiary amines (no H) and quaternary ammonium
    monoamine_pattern = Chem.MolFromSmarts("[a;r]:[#6](-[#6]-[N&!H0&!$(N-[-0])])")  # Two-carbon chain, amino with H
    matches = mol.GetSubstructMatches(monoamine_pattern)

    # Check for exactly one amino group meeting criteria
    if len(matches) != 1:
        return False, f"Found {len(matches)} amino groups with two-carbon chain to aromatic ring"

    # Verify the connecting chain is exactly two carbons
    for match in matches:
        aromatic_idx, carbon1, carbon2, nitrogen_idx = match
        # Ensure the two carbons are in a chain without branching
        atom1 = mol.GetAtomWithIdx(carbon1)
        atom2 = mol.GetAtomWithIdx(carbon2)
        if atom1.GetDegree() != 2 or atom2.GetDegree() != 2:
            return False, "Branched chain detected"

    return True, "Amino group connected to aromatic ring via two-carbon chain"