"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound involves the presence of iodine covalently bonded to carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification or misclassification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    has_iodine_carbon_bond = False

    # Traverse through the atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 53:  # Atomic number for iodine
            # Check neighboring atoms to find carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Atomic number for carbon
                    has_iodine_carbon_bond = True
                    break

    if has_iodine_carbon_bond:
        return True, "Contains at least one iodine-carbon bond, typical of organoiodine compounds."

    return False, "No iodine-carbon bonds found, not an organoiodine compound."