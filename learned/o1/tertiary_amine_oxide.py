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
    A tertiary amine oxide has a nitrogen atom with a formal +1 charge,
    bonded to three organic groups (typically carbon atoms) and an oxygen atom with a -1 charge.

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
        # Check if the atom is nitrogen
        if atom.GetAtomicNum() == 7:
            # Check if nitrogen has formal charge +1
            if atom.GetFormalCharge() != +1:
                continue
            # Check if nitrogen has zero hydrogens (no N-H bonds)
            if atom.GetTotalNumHs() != 0:
                continue
            # Check if nitrogen is connected to exactly 4 atoms (3 groups + 1 oxygen)
            if atom.GetDegree() != 4:
                continue

            # Initialize counters
            num_carbons = 0
            num_oxygens = 0
            has_oxygen_with_neg_charge = False

            # Iterate over neighbors of nitrogen
            for nbr in atom.GetNeighbors():
                atomic_num = nbr.GetAtomicNum()
                # Count carbon neighbors
                if atomic_num == 6:
                    num_carbons += 1
                # Check for oxygen neighbor
                elif atomic_num == 8:
                    num_oxygens += 1
                    # Check if oxygen has formal charge -1
                    if nbr.GetFormalCharge() == -1:
                        has_oxygen_with_neg_charge = True
                else:
                    # Nitrogen is bonded to an atom other than carbon or oxygen
                    break
            else:
                # Check if nitrogen is connected to exactly 3 carbons and 1 oxygen with -1 charge
                if num_carbons == 3 and num_oxygens == 1 and has_oxygen_with_neg_charge:
                    return True, "Contains tertiary amine oxide functional group"
    return False, "Does not contain tertiary amine oxide functional group"

__metadata__ = {   'chemical_class': {   'name': 'tertiary amine oxide',
                                         'definition': 'An N-oxide where there are three organic groups bonded to the nitrogen atom.'}}