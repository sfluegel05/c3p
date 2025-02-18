"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: CHEBI:XXXXX hydroxynaphthoquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone with at least one hydroxy group attached to the naphthalene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general naphthoquinone pattern (matches both 1,2 and 1,4-naphthoquinones)
    naphthoquinone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)C(=O)C(=O)C=C2")
    if not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No naphthoquinone moiety found"

    # Look for at least one hydroxy group attached to the naphthalene ring
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) == 0:
        return False, "No hydroxy group found"

    # Ensure the hydroxy group is directly attached to the naphthalene ring
    naphthalene_atoms = set()
    for match in mol.GetSubstructMatches(naphthoquinone_pattern):
        naphthalene_atoms.update(match)

    for hydroxy_match in hydroxy_matches:
        hydroxy_atom = mol.GetAtomWithIdx(hydroxy_match[0])
        neighbor_atom = hydroxy_atom.GetNeighbors()[0]
        if neighbor_atom.GetIdx() not in naphthalene_atoms:
            return False, "Hydroxy group not directly attached to the naphthalene ring"

    return True, "Contains a naphthoquinone moiety with at least one hydroxy group directly attached to the naphthalene ring"