"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is a nitrogen atom with a positive charge and two directly attached
    carbon containing groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1: # Check for nitrogen and positive charge
            carbon_neighbors = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6: #Check if neighbour is a carbon
                    carbon_neighbors += 1
            if carbon_neighbors == 2:
                return True, "Contains at least one positively charged nitrogen with two directly attached carbon groups."


    return False, "Does not contain a positively charged nitrogen with two directly attached carbon groups."