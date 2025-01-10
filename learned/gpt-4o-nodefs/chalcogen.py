"""
Classifies: CHEBI:33303 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a chemical entity is a chalcogen based on its SMILES string.
    Chalcogens include atoms of oxygen (O), sulfur (S), selenium (Se), tellurium (Te), and polonium (Po),
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

    # Check each atom in the molecule
    for atom in mol.GetAtoms():
        # Extract the chemical symbol and isotope if any
        symbol = atom.GetSymbol()
        isotope = atom.GetIsotope()
        
        # Concrete isotope check (e.g. [16O]) will have atomic number matching
        if symbol in chalcogen_elements:
            return True, f"Chemical entity is a chalcogen: {symbol} isotope {isotope}"

    return False, "Chemical entity is not a chalcogen"