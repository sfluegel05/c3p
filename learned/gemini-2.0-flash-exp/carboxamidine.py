"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine has the structure RC(=NR)NR2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, check if it is guanidine C(=N)(N)N, and exclude.
    guanidine_pattern = Chem.MolFromSmarts("[C](=[N])([N])-[N]")
    if mol.HasSubstructMatch(guanidine_pattern):
         return False, "Molecule is a Guanidine"

    # Check for core carboxamidine pattern [C](=[N])[N]
    # The [C] must be doubly bonded to [N], the other [N] must be singly bonded.
    # Avoid Nitro groups (N(=O)=O), nitriles (C#N) and azides (N=N=N)
    carboxamidine_pattern = Chem.MolFromSmarts("[CX3]=[NX1]-[NX2]")
    matches = mol.GetSubstructMatches(carboxamidine_pattern)

    if not matches:
        return False, "Core carboxamidine structure C(=N)-N not found."

    return True, "Molecule matches the core structure of a carboxamidine"