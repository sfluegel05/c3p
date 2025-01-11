"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion arises from the protonation of a secondary amine, displaying
    a positively charged nitrogen bonded to exactly two carbon atoms and a hydrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if the atom is a positively charged nitrogen (ammonium ion)
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Get the neighbors of the nitrogen atom
            neighbors = atom.GetNeighbors()
            carbon_count = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() == 6)

            # Check for secondary amine configuration: 2 carbons
            if carbon_count == 2:
                return True, "Contains protonated secondary amine group forming secondary ammonium ion"

    return False, "Does not contain the features of a secondary ammonium ion"