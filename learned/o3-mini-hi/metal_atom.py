"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: A metal atom.
A metal atom is defined as an atom that exhibits typical metallic properties,
being typically shiny with high electrical and thermal conductivity.
This function determines whether a given SMILES string represents a neutral metal atom.
The classification requires that the molecule has exactly one atom and that the atom is neutral.
"""

from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    The classification criteria are:
      - The molecule must be a single atom.
      - The atom must have a formal charge of 0 (i.e. must be neutral).
      - The atomic symbol must be in the list of known metal elements.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a neutral metal atom, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # The molecule should have exactly one atom
    if mol.GetNumAtoms() != 1:
        return False, "The molecule contains more than one atom; expected a single atom."

    # Retrieve the single atom
    atom = mol.GetAtomWithIdx(0)
    
    # Check that the atom is neutral (formal charge must be 0)
    formal_charge = atom.GetFormalCharge()
    if formal_charge != 0:
        return False, f"The atom has a formal charge of {formal_charge}, expected a neutral atom."

    # Retrieve the atomic symbol (ignoring isotope information)
    symbol = atom.GetSymbol()

    # Define a set of metal elements (includes alkali metals, alkaline earth metals,
    # transition metals, post-transition metals, lanthanides, actinides, and superheavy metals).
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

    # Check if the atom's symbol is recognized as a metal
    if symbol in metal_elements:
        return True, f"The atom is a metal: {symbol}"
    else:
        return False, f"Atom symbol '{symbol}' is not recognized as a metal"

# Example usage (uncomment for testing):
# print(is_metal_atom("[No]"))      # Expected: True (nobelium, neutral)
# print(is_metal_atom("[Mn+7]"))     # Expected: False (charged manganese species)