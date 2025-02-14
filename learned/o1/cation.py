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
    
    # Iterate over all atoms to check for positive formal charges
    has_positive_charge = False
    for atom in mol.GetAtoms():
        formal_charge = atom.GetFormalCharge()
        if formal_charge > 0:
            has_positive_charge = True
            break  # No need to check further if a positive charge is found

    if has_positive_charge:
        return True, "Molecule has one or more positive charges; it is a cation"
    else:
        return False, "Molecule has no positive charges; it is not a cation"