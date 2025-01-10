"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Attempts to determine if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    Dihydroagarofuran sesquiterpenoids are complex structures with specific ring systems and functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined pattern to match more accurately typical dihydroagarofuran core features
    # This example assumes a typical bicyclic framework for dihydroagarofuran cores
    dihydroagarofuran_core_pattern = Chem.MolFromSmarts("C1[C@@H]2[C@H](OC(=O)C)[C@@H]3OC(=O)[C@H]23")  # Example refined pattern

    if not mol.HasSubstructMatch(dihydroagarofuran_core_pattern):
        return False, "Core structure of dihydroagarofuran not found"

    # Check for more specific functional group / ester presence
    # Typical dihydroagarofuran sesquiterpenoids have extensive ester groups
    ester_pattern = Chem.MolFromSmarts("O=C(O)C")  # Basic ester pattern
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:  # Expect at least two ester matches in these molecules
        return False, f"Ester groups absent or fewer than expected, found {len(ester_matches)}"

    return True, "Molecule matches dihydroagarofuran sesquiterpenoid core structure with ester characteristics"

# Test examples (update patterns based on dihydroagarofuran specific knowledge if necessary)