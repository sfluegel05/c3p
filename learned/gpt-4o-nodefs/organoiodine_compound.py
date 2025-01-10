"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound is characterized by iodine covalently bonded to carbon,
    with specific structural arrangements typical of this class.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification or misclassification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    has_iodine_carbon_bond = False
    iodine_count = 0
    # Traverse through the atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 53:  # Atomic number for iodine
            iodine_count += 1
            # Check neighboring atoms to find carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Atomic number for carbon
                    has_iodine_carbon_bond = True
                    break

    # Check if key structural criteria for organoiodine compounds are met
    if has_iodine_carbon_bond and iodine_count > 1:
        return True, "Iodine is bonded directly to carbon, indicating an organoiodine compound with multiple Iodine atoms."

    return False, "No characteristic structure of organoiodine compound found."