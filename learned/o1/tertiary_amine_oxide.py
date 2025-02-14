"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: tertiary amine oxide
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide has a nitrogen atom bonded to an oxygen atom,
    and the nitrogen is also bonded to three organic groups (typically carbon atoms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Check if nitrogen is connected to exactly 4 atoms
            if atom.GetDegree() != 4:
                continue
            neighbors = atom.GetNeighbors()
            num_carbons = 0
            num_oxygens = 0
            other_atoms = 0
            for nbr in neighbors:
                atomic_num = nbr.GetAtomicNum()
                if atomic_num == 6:
                    num_carbons += 1
                elif atomic_num == 8:
                    num_oxygens += 1
                else:
                    other_atoms +=1
            # Check if nitrogen is connected to exactly 3 carbons and 1 oxygen
            if num_carbons == 3 and num_oxygens == 1 and other_atoms == 0:
                return True, "Contains tertiary amine oxide functional group"
    return False, "Does not contain tertiary amine oxide functional group"

__metadata__ = {   'chemical_class': {   'name': 'tertiary amine oxide',
                                         'definition': 'An N-oxide where there are three organic groups bonded to the nitrogen atom.'}}