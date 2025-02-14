"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is a monoatomic or polyatomic species having one or more elementary charges of the proton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cation, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if any atom has a positive formal charge
    has_positive_charge = any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms())

    if has_positive_charge:
        return True, "Contains atoms with positive formal charges; it is a cation"
    else:
        return False, "No atoms with positive formal charges; it is not a cation"