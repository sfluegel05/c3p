"""
Classifies: CHEBI:33313 polonium atom
"""
"""
Classifies: CHEBI:33313 polonium atom
"""
from rdkit import Chem

def is_polonium_atom(smiles: str):
    """
    Determines if a molecule is a polonium atom based on its SMILES string.
    A polonium atom is represented as [Po] or [<isotope>Po] in SMILES.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polonium atom, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule consists of a single atom
    if mol.GetNumAtoms() != 1:
        return False, "Molecule must consist of a single atom"

    # Get the single atom
    atom = mol.GetAtomWithIdx(0)

    # Check if the atom is polonium (atomic number 84)
    if atom.GetAtomicNum() != 84:
        return False, "Atom is not polonium"

    # Check if the atom has an isotope specification
    isotope = atom.GetIsotope()
    if isotope != 0:
        # Polonium isotopes range from 188 to 220
        if isotope < 188 or isotope > 220:
            return False, f"Invalid polonium isotope: {isotope}"

    return True, "Single polonium atom, possibly with an isotope specification"