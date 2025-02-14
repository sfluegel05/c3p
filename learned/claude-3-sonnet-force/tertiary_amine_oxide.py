"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: CHEBI:50786 tertiary amine oxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide is an N-oxide where there are three organic groups bonded to the nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N-oxide patterns ([N+](=O)[!O;!H] or [N+](=[O-])[C])
    n_oxide_patterns = [Chem.MolFromSmarts("[N+](=O)[!O;!H]"), Chem.MolFromSmarts("[N+](=[O-])[C]")]
    n_oxide_matches = []
    for pattern in n_oxide_patterns:
        n_oxide_matches.extend(mol.GetSubstructMatches(pattern))
    if len(n_oxide_matches) != 1:
        return False, f"Found {len(n_oxide_matches)} N-oxide groups, need exactly one"

    # Check if N-oxide is tertiary (3 organic groups attached)
    n_oxide_idx = n_oxide_matches[0][0]
    n_oxide_atom = mol.GetAtomWithIdx(n_oxide_idx)
    organic_neighbors = [nbr for nbr in n_oxide_atom.GetNeighbors() if nbr.GetAtomicNum() != 8 and nbr.GetAtomicNum() != 1]
    if len(organic_neighbors) != 3:
        return False, f"N-oxide has {len(organic_neighbors)} organic neighbors, need exactly 3"

    # Check if organic groups are alkyl/aryl
    alkyl_aryl_pattern = Chem.MolFromSmarts("[C;A]")
    for nbr in organic_neighbors:
        if not mol.GetAtomWithIdx(nbr.GetIdx()).HasSubstructMatch(alkyl_aryl_pattern):
            return False, "One of the organic groups is not alkyl or aryl"

    # Check if N-oxide is part of a ring system
    ring_info = mol.GetRingInfo()
    is_ring_atom = ring_info.IsAtomRingBond(n_oxide_idx)
    if is_ring_atom:
        # Implement additional checks for ring systems, e.g., ring size, aromaticity, etc.
        pass

    # Check for additional substituents/functional groups on organic groups
    for nbr in organic_neighbors:
        nbr_atom = mol.GetAtomWithIdx(nbr.GetIdx())
        for atom in nbr_atom.GetNeighbors():
            if atom.GetAtomicNum() not in [1, 6, 7, 8, 9, 15, 16, 17, 35]:
                # Implement additional checks for specific substituents/functional groups
                pass

    # Check stereochemistry and connectivity
    # ...

    return True, "Molecule is a tertiary amine oxide"