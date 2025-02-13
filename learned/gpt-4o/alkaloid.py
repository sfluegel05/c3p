"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determine if a molecule is likely an alkaloid based on its SMILES string.
    This function checks for basic nitrogen in a heterocyclic structure and ensures
    it is not merely an exocyclic nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nitrogen atoms in the molecule
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atom found"

    # Check if any nitrogen is part of a heterocyclic ring
    is_heterocyclic = False
    for atom in nitrogen_atoms:
        if atom.IsInRing():
            is_heterocyclic = True
            break

    if not is_heterocyclic:
        return False, "No nitrogen atom found in a ring structure (heterocyclic)"

    # Ensure nitrogen is not solely exocyclic
    exocyclic_nitrogens = [atom for atom in nitrogen_atoms if not atom.IsInRing()]
    if len(exocyclic_nitrogens) == len(nitrogen_atoms):
        return False, "All nitrogen atoms are exocyclic"

    # If molecular features fit the broad definition, classify as alkaloid
    return True, "Contains basic nitrogen in heterocyclic structure, likely an alkaloid"

# Example usage:
# smiles = "CCCCCCCCCCC1=CC=C(N1)\C=C1/N=C(C=C1OC)C1=CC=CN1"  # Example SMILES
# result, reason = is_alkaloid(smiles)
# print(result, reason)