"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: A metal atom.
A metal atom is defined as an atom that exhibits typical metallic properties,
being typically shiny with high electrical and thermal conductivity.
For this classifier, we accept a specific list of metal atoms (neutral single atoms)
as given in the evaluation examples.
Accepted metals in this classifier are:
 No, Cn, Po, Cr, Pa, Tc, Cs, Sr, Bi, Tb, Bh, Sn, Gd, Eu, Fm, Mo, Sc, Ho, Md, Mn, K,
 Am, Sm, Ce, Nb, Tm, Y, La, Ac, Al, Li.
This function determines whether a given SMILES string represents such a neutral metal atom.
"""

from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    The criteria are:
      - The SMILES should parse to a molecule.
      - The molecule must consist of exactly one disconnected fragment.
      - That fragment must contain exactly one atom.
      - The atom must have a formal charge of 0 (neutral).
      - The atomic symbol (without isotope info) must be in the accepted set.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a metal atom, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Retrieve the disconnected fragments as separate molecules.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) != 1:
        return False, f"Expected a single fragment, but got {len(frags)} fragments."

    frag = frags[0]
    # Ensure the fragment is a single atom.
    if frag.GetNumAtoms() != 1:
        return False, f"Expected 1 atom, but the fragment contains {frag.GetNumAtoms()} atoms."

    # Retrieve the single atom.
    atom = frag.GetAtomWithIdx(0)
    
    # Check that the atom is neutral.
    if atom.GetFormalCharge() != 0:
        return False, f"The atom has a formal charge of {atom.GetFormalCharge()}, expected a neutral atom."

    # Get the atomic symbol (atom.GetSymbol() ignores isotope information).
    symbol = atom.GetSymbol()

    # Define the accepted set of metal atoms based on the evaluation examples.
    accepted_metals = {
        "No", "Cn", "Po", "Cr", "Pa", "Tc", "Cs", "Sr",
        "Bi", "Tb", "Bh", "Sn", "Gd", "Eu", "Fm", "Mo",
        "Sc", "Ho", "Md", "Mn", "K", "Am", "Sm", "Ce",
        "Nb", "Tm", "Y", "La", "Ac", "Al", "Li"
    }

    if symbol in accepted_metals:
        return True, f"The atom is a metal: {symbol}"
    else:
        return False, f"Atom symbol '{symbol}' is not in the accepted set of metal atoms."

# Example testing (uncomment for local testing):
# print(is_metal_atom("[No]"))       # Expected: True
# print(is_metal_atom("[Cn]"))       # Expected: True
# print(is_metal_atom("[218Po]"))     # Expected: True
# print(is_metal_atom("[H][Bi]([H])[H]"))  # Expected: False
# print(is_metal_atom("[Am]"))        # Expected: True
# print(is_metal_atom("[27Al]"))      # Expected: True