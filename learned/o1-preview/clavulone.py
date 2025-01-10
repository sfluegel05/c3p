"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    A clavulone is an esterified prostanoid obtained from marine corals, generally characterized by a substituted cyclopentenone ring, conjugated systems, ester groups, long aliphatic chains, and often halogen substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for substituted cyclopentenone ring with conjugation
    cyclopentenone_pattern = Chem.MolFromSmarts('O=C1C=CC=C1')  # Conjugated cyclopentenone
    if not mol.HasSubstructMatch(cyclopentenone_pattern):
        return False, "No conjugated cyclopentenone ring found"

    # Check for ester groups (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts('[$([CX3](=O)[OX2H0])]')  # Ester group
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found"

    # Check for long aliphatic chains (chain of at least 6 carbons)
    chain_pattern = Chem.MolFromSmarts('C' + ('C' * 5) + '[C,C+]')  # Chain of at least 6 carbons
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long aliphatic chain of at least 6 carbons found"

    # Check for halogen substituents (Cl, Br, I) attached to the ring
    halogen_pattern = Chem.MolFromSmarts('[#6]-[Cl,Br,I]')
    if not mol.HasSubstructMatch(halogen_pattern):
        return False, "No halogen substituents found attached to carbon"

    # Ensure that all key features are connected appropriately
    # Combine cyclopentenone, ester, and halogen patterns
    combined_pattern = Chem.MolFromSmarts('[$(O=C1C=CC=C1)]$([#6]-[Cl,Br,I])[$([CX3](=O)[OX2H0])])')
    if not mol.HasSubstructMatch(combined_pattern):
        return False, "Key features are not connected appropriately"

    # Check if molecule has at least one ring (prostanoid core)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 1:
        return False, "Too few rings for a clavulone"

    return True, "Molecule matches features characteristic of clavulones"

__metadata__ = {
    'chemical_class': {
        'name': 'clavulone',
        'definition': 'A class of esterified prostanoids obtained from marine corals.'
    }
}