"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: organohalogen compound
Definition: 'A compound containing at least one carbon-halogen bond (where X is a halogen atom).'
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound contains at least one carbon-halogen bond (C-X), where X is F, Cl, Br, or I.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Atomic numbers for carbon and halogens
    carbon_atomic_num = 6
    halogen_atomic_nums = [9, 17, 35, 53]  # F, Cl, Br, I

    # Iterate over all bonds in the molecule
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atom1_num = atom1.GetAtomicNum()
        atom2_num = atom2.GetAtomicNum()

        # Check if one atom is carbon and the other is a halogen
        if ((atom1_num == carbon_atomic_num and atom2_num in halogen_atomic_nums) or
            (atom2_num == carbon_atomic_num and atom1_num in halogen_atomic_nums)):
            return True, "Contains at least one carbon-halogen bond"

    # No carbon-halogen bonds found
    return False, "No carbon-halogen bonds found"