"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: A metal atom.
A metal atom is defined as an atom of an element that exhibits typical metallic properties,
being typically shiny, with high electrical and thermal conductivity.
This function determines whether a given SMILES string represents a metal atom.
"""

from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    A valid metal atom SMILES is a single-atom molecule and that atom is one of the elements
    known to display metallic properties.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a metal atom, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that the molecule has exactly one atom
    if mol.GetNumAtoms() != 1:
        return False, "The molecule contains more than one atom; expected a single atom."

    # Retrieve the atom symbol (this ignores isotope information; e.g., [218Po] will yield 'Po')
    atom = mol.GetAtomWithIdx(0)
    symbol = atom.GetSymbol()

    # Define a set of metal elements.
    # This list includes the alkali metals, alkaline earth metals, most transition metals,
    # post-transition metals, lanthanides, and actinides.
    metal_elements = {
        # Alkali metals
        "Li", "Na", "K", "Rb", "Cs", "Fr",
        # Alkaline earth metals
        "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
        # Transition metals
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        # Post-transition metals
        "Al", "Ga", "In", "Sn", "Tl", "Pb", "Bi", "Po",
        # Lanthanides
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
        "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        # Actinides
        "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
        "Es", "Fm", "Md", "No", "Lr",
        # Recently named or superheavy elements considered metals
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
        "Fl", "Mc", "Lv", "Ts", "Og"
    }

    # Check if the atom's symbol is in the set of metal elements.
    if symbol in metal_elements:
        return True, f"The atom is a metal: {symbol}"
    else:
        return False, f"Atom symbol '{symbol}' is not recognized as a metal"

# Example usage (uncomment for testing):
# print(is_metal_atom("[No]"))       # Expected: True (nobelium)
# print(is_metal_atom("[C]"))         # Expected: False (carbon is nonmetal)
# print(is_metal_atom("CC"))          # Expected: False (more than one atom)