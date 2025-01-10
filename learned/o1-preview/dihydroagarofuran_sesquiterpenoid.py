"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: dihydroagarofuran sesquiterpenoid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid is a sesquiterpenoid (15-carbon terpenoid) with a dihydroagarofuran skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the dihydroagarofuran skeleton as a SMARTS pattern
    # The skeleton consists of a tricyclic framework with a tetrahydrofuran ring fused to two cyclopentane rings
    dihydroagarofuran_smarts = """
    [C@@H]1[C@H]2[C@@H](C[C@H](O)C1)C[C@@H]3[C@@H](C2)CC3
    """
    skeleton = Chem.MolFromSmarts(dihydroagarofuran_smarts)
    if skeleton is None:
        return False, "Invalid SMARTS pattern for dihydroagarofuran skeleton"

    # Check if the molecule contains the dihydroagarofuran skeleton
    if not mol.HasSubstructMatch(skeleton):
        return False, "Dihydroagarofuran skeleton not found"

    # Check for sesquiterpenoid (15-carbon terpenoid)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 15:
        return False, f"Number of carbons ({num_carbons}) less than 15"

    # Optional: Check for terpenoid-like oxygenation patterns or substituents if necessary
    # For example, check for ester or hydroxyl groups common in sesquiterpenoids
    num_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygen < 1:
        return False, "No oxygen atoms found, unlikely to be a dihydroagarofuran sesquiterpenoid"

    return True, "Contains dihydroagarofuran skeleton and meets criteria for dihydroagarofuran sesquiterpenoid"