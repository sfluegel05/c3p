"""
Classifies: CHEBI:33303 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a chemical entity is a chalcogen based on its SMILES string.
    Chalcogens include standalone atoms of oxygen (O), sulfur (S), selenium (Se), tellurium (Te), and polonium (Po),
    as well as their isotopes.

    Args:
        smiles (str): SMILES string of the chemical entity

    Returns:
        bool: True if the chemical entity is a chalcogen, False otherwise
        str: Reason for classification
    """

    # Define chalcogen elements
    chalcogen_elements = {"O", "S", "Se", "Te", "Po"}

    # Parse SMILES to obtain the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule consists of exactly one atom
    if mol.GetNumAtoms() == 1:
        atom = mol.GetAtomWithIdx(0)
        symbol = atom.GetSymbol()
        if symbol in chalcogen_elements:
            isotope = atom.GetIsotope()
            return True, f"Chemical entity is a chalcogen: {symbol} isotope {isotope}"

    return False, "Chemical entity is not a chalcogen"