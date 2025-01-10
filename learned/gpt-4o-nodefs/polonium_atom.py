"""
Classifies: CHEBI:33313 polonium atom
"""
from rdkit import Chem

def is_polonium_atom(smiles: str):
    """
    Determines if a molecule is a polonium atom based on its SMILES string.
    A polonium atom in SMILES is represented as [massPo], where 'mass' represents
    its atomic mass number.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polonium atom, False otherwise
        str: Reason for classification
    """

    # Parse SMILES and catch errors
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
    except:
        return False, "Error parsing SMILES"

    # Look for polonium atom with isotopic specification
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'Po':
            isotope = atom.GetIsotope()
            if isotope != 0:  # Isotope should be non-zero to be specified
                return True, f"Contains polonium isotope: [{isotope}Po]"
            else:
                return False, "Polonium found without isotope specification"

    return False, "No polonium atom found in SMILES"