"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: CHEBI:38194 cannabinoid
A diverse group of pharmacologically active secondary metabolite characteristic to Cannabis plant
as well as produced naturally in the body by humans and animals. Cannabinoids contain oxygen as a
part of the heterocyclic ring or in the form of various functional groups. They are subdivided on
the basis of their origin.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_cannabinoid(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of a pyran or benzopyran ring
    # Cannabinoids typically contain a pyran or benzopyran ring
    pyran_pattern = Chem.MolFromSmarts("O1CCCCC1")
    benzopyran_pattern = Chem.MolFromSmarts("c1ccoc2ccccc12")
    if not mol.HasSubstructMatch(pyran_pattern) and not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No pyran or benzopyran ring found"

    # Check for long aliphatic chains
    # Cannabinoids typically have long aliphatic chains
    alkyl_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC")
    alkyl_chain_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    if not alkyl_chain_matches:
        return False, "No long aliphatic chain found"

    # Check for presence of oxygen-containing functional groups
    # Cannabinoids contain oxygen in various functional groups
    oxygen_functional_groups = ['O', 'C(=O)', 'C(=O)O', 'OC(=O)']
    oxygen_functional_group_pattern = Chem.MolFromSmarts("|".join(oxygen_functional_groups))
    if not mol.HasSubstructMatch(oxygen_functional_group_pattern):
        return False, "No oxygen-containing functional group found"

    return True, "Contains a pyran or benzopyran ring, long aliphatic chains, and oxygen-containing functional groups, which are characteristic of cannabinoids"