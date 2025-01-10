"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    A clavulone is an esterified prostanoid obtained from marine corals, generally characterized by a cyclopentenone ring with conjugated double bonds, ester groups, long aliphatic chains, and often halogen substitutions.

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

    # Check for cyclopentenone ring (five-membered ring with ketone and conjugated double bonds)
    cyclopentenone_pattern = Chem.MolFromSmarts('O=C1C=CC=C1')  # Simplified pattern for cyclopentenone
    if not mol.HasSubstructMatch(cyclopentenone_pattern):
        return False, "No cyclopentenone ring found"

    # Check for ester groups (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts('C(=O)O')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester groups found"

    # Check for long aliphatic chains attached via conjugated double bonds
    # Simplified as chain of at least 6 carbons with double bonds
    long_chain_pattern = Chem.MolFromSmarts('C=CC=CC=CC')  # Pattern for conjugated chain
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long conjugated aliphatic chain found"

    # Check for halogen substituents on the ring (Cl, Br, I)
    halogen_pattern = Chem.MolFromSmarts('[#6;r5][F,Cl,Br,I]')
    if not mol.HasSubstructMatch(halogen_pattern):
        return False, "No halogen substituent on the ring found"

    # Check if molecule is an esterified prostanoid
    # Prostanoids generally have specific molecular features
    # For simplification, we can check for molecular weight and number of rings
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a clavulone"

    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 1:
        return False, "Too few rings for a clavulone"

    return True, "Molecule matches features characteristic of clavulones"

__metadata__ = {   'chemical_class': {   'name': 'clavulone',
                              'definition': 'A class of esterified prostanoids obtained from marine corals.'}
}