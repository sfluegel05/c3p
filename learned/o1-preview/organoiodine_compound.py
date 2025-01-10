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
    An organoiodine compound contains at least one carbon-iodine (C-I) bond where the carbon is
    bonded to iodine and at least one other non-halogen atom (excluding hydrogen).

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

    # Iterate over all bonds
    for bond in mol.GetBonds():
        # Check if bond is between carbon and iodine
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if ((atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 53) or
            (atom1.GetAtomicNum() == 53 and atom2.GetAtomicNum() == 6)):
            # Identify the carbon atom
            carbon_atom = atom1 if atom1.GetAtomicNum() == 6 else atom2
            iodine_atom = atom1 if atom1.GetAtomicNum() == 53 else atom2
            # Check if carbon is bonded to at least one non-halogen atom (excluding hydrogen)
            non_halogen_neighbors = False
            for neighbor in carbon_atom.GetNeighbors():
                neigh_atomic_num = neighbor.GetAtomicNum()
                if neighbor.GetIdx() != iodine_atom.GetIdx():
                    if neigh_atomic_num not in [1, 9, 17, 35, 53]:
                        non_halogen_neighbors = True
                        break
            if non_halogen_neighbors:
                return True, "Contains carbon-iodine bond where carbon is bonded to at least one non-halogen atom"
    return False, "No suitable carbon-iodine bonds found"