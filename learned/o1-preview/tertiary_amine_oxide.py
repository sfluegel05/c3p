"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: tertiary amine oxide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide is an N-oxide where the nitrogen atom is bonded to three organic groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize flags
    is_tertiary_amine_oxide = False

    # Iterate over atoms to find positively charged nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == +1 and atom.GetDegree() == 4:
            # Initialize counters
            oxygen_found = False
            non_hydrogen_neighbors = []

            # Iterate over neighbors of nitrogen atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetFormalCharge() == -1:
                    # Found oxygen with negative charge (oxide)
                    oxygen_found = True
                elif neighbor.GetAtomicNum() == 1:
                    # Neighbor is hydrogen, not an organic group
                    break  # Not a tertiary amine oxide, move to next nitrogen atom
                else:
                    # Add neighbor to list of non-hydrogen neighbors
                    non_hydrogen_neighbors.append(neighbor)
            else:
                # Ensure exactly one oxygen and three other neighbors
                if oxygen_found and len(non_hydrogen_neighbors) == 3:
                    is_tertiary_amine_oxide = True
                    break  # Found a tertiary amine oxide

    if is_tertiary_amine_oxide:
        return True, "Contains a positively charged nitrogen bonded to an oxygen (-1 charge) and three organic groups"
    else:
        return False, "Does not contain the tertiary amine oxide functional group with three organic groups attached to nitrogen"

__metadata__ = {   'chemical_class': {   'name': 'tertiary amine oxide',
                                         'definition': 'An N-oxide where there are three organic groups bonded to the nitrogen atom.'}}