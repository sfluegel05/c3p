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
    # Pattern: aromatic atom -> any carbon -> any carbon -> nitrogen (any substitution except nitro/amide)
    # Excludes nitro groups and amides by not matching those patterns
    monoamine_pattern = Chem.MolFromSmarts("[a]~[#6]~[#6]~[N]")
    matches = mol.GetSubstructMatches(monoamine_pattern)

    # Check for exactly one amino group meeting criteria
    if len(matches) != 1:
        return False, f"Found {len(matches)} amino groups with two-carbon chain to aromatic ring"

    # Additional check to exclude nitro groups and amides
    # Nitro pattern: [N+](=O)[O-]
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    if mol.HasSubstructMatch(nitro_pattern):
        return False, "Nitro group present"

    # Amide pattern: [NX3][CX3](=O)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Amide group present"

    return True, "Amino group connected to aromatic ring via two-carbon chain"