"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:36866 carotenoid
One of a class of tetraterpenoids (C40), formally derived from the acyclic parent,
psi,psi-carotene by hydrogenation, dehydrogenation, cyclization, oxidation, or
combination of these processes. This class includes carotenes, xanthophylls and
certain compounds that arise from rearrangement of the skeleton of psi,psi-carotene
or by loss of part of this structure. Retinoids are excluded.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for C40 backbone
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 40:
        return False, "Not a C40 compound"

    # Check for tetraterpenoid scaffold
    terpene_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]([C@@H]([C@@H]([C@H](C1)C)C)C)C")
    terpene_matches = mol.GetSubstructMatches(terpene_pattern)
    if len(terpene_matches) < 4:
        return False, "Not a tetraterpenoid scaffold"

    # Check for cyclization, dehydrogenation, oxidation patterns
    cyclized_pattern = Chem.MolFromSmarts("C1CCCCC1")
    dehydrogenated_pattern = Chem.MolFromSmarts("C=C")
    oxidized_pattern = Chem.MolFromSmarts("O")

    cyclized_matches = mol.GetSubstructMatches(cyclized_pattern)
    dehydrogenated_matches = mol.GetSubstructMatches(dehydrogenated_pattern)
    oxidized_matches = mol.GetSubstructMatches(oxidized_pattern)

    if not (cyclized_matches or dehydrogenated_matches or oxidized_matches):
        return False, "No evidence of cyclization, dehydrogenation or oxidation"

    # Check for absence of retinoid scaffold
    retinoid_pattern = Chem.MolFromSmarts("CC(C)(C)c1ccc(O)cc1")
    retinoid_matches = mol.GetSubstructMatches(retinoid_pattern)
    if retinoid_matches:
        return False, "Contains retinoid scaffold"

    return True, "Molecule is a carotenoid"