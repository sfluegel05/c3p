"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: organoiodine compound
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound is defined as a compound containing at least one carbon-iodine (C-I) bond,
    where the carbon is bonded to at least one carbon or heteroatom other than halogens (F, Cl, Br, I) and hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over atoms to find carbon-iodine bonds
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        # Check if the bond is between carbon and iodine
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 53) or \
           (atom1.GetAtomicNum() == 53 and atom2.GetAtomicNum() == 6):

            # Identify the carbon atom
            if atom1.GetAtomicNum() == 6:
                carbon = atom1
            else:
                carbon = atom2

            # Check neighbors of the carbon atom
            neighbor_atomic_nums = [neighbor.GetAtomicNum() for neighbor in carbon.GetNeighbors()]

            # Exclude halogens (F=9, Cl=17, Br=35, I=53) and hydrogen (H=1)
            non_halogens = [num for num in neighbor_atomic_nums if num not in [1, 9, 17, 35, 53]]

            if len(non_halogens) > 0:
                return True, "Contains carbon-iodine bond where carbon is bonded to non-halogen atom"

    return False, "No suitable carbon-iodine bonds found or carbon only bonded to halogens/hydrogen"