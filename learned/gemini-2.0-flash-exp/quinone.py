"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is a cyclic dione structure derived from an aromatic compound, with conjugated double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a quinone substructure WITHIN a ring
    # This pattern looks for a 6-membered aromatic ring with at least two carbonyl groups (=O)
    # separated by a conjugated system. We are also explicitly requiring that the atoms are aromatic (lowercase c).
    # The two carbonyl carbons are explicitly required to be within the same ring.
    quinone_pattern = Chem.MolFromSmarts("[c:1]1[c:2][c:3][c:4](=[O])[c:5][c:6](=[O])1")
    quinone_pattern_5ring = Chem.MolFromSmarts("[c:1]1[c:2][c:3](=[O])[c:4][c:5](=[O])1")

    # Check for the presence of the quinone pattern
    if mol.HasSubstructMatch(quinone_pattern) or mol.HasSubstructMatch(quinone_pattern_5ring):
        return True, "Aromatic ring with dione substructure detected"

    return False, "No quinone structure detected"