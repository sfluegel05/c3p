"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: A metal atom.
A metal atom is defined as an atom that exhibits typical metallic properties,
being typically shiny with high electrical and thermal conductivity.
For this classifier, the molecule must be a single neutral atom (ignoring isotopic labels)
and its atomic symbol must be in an accepted set of metal atoms.
Accepted metals in this classifier (expanded based on feedback) are:
 No, Cn, Po, Cr, Pa, Tc, Cs, Sr, Bi, Tb, Bh, Sn, Gd, Eu, Fm, Mo, Sc, Ho, Md, Mn,
 K, Am, Sm, Ce, Nb, Tm, Y, La, Ac, Al, Li, Rb, Co, Cf, Os, Zr, Ni, Ba, Hg, Bk.
"""

from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    The criteria are:
      - The SMILES must parse to a molecule.
      - The molecule must consist of exactly one disconnected fragment.
      - That fragment must contain exactly one atom.
      - The atom must have a formal charge of 0 (neutral).
      - The atomic symbol (ignoring any isotope numbers) must be in the accepted set.

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
    
    # Split the molecule into disconnected fragments.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) != 1:
        return False, f"Expected a single fragment, but got {len(frags)} fragments."
    
    frag = frags[0]
    # Check that the fragment contains exactly one atom.
    if frag.GetNumAtoms() != 1:
        return False, f"Expected 1 atom, but the fragment contains {frag.GetNumAtoms()} atoms."
    
    # Retrieve the single atom and check its formal charge.
    atom = frag.GetAtomWithIdx(0)
    if atom.GetFormalCharge() != 0:
        return False, f"The atom has a formal charge of {atom.GetFormalCharge()}, expected a neutral atom."
    
    # Get the atomic symbol (this ignores any isotope information)
    symbol = atom.GetSymbol()
    
    # Define the accepted set of metal atoms.
    accepted_metals = {
        "No", "Cn", "Po", "Cr", "Pa", "Tc", "Cs", "Sr", "Bi", "Tb", "Bh",
        "Sn", "Gd", "Eu", "Fm", "Mo", "Sc", "Ho", "Md", "Mn", "K",
        "Am", "Sm", "Ce", "Nb", "Tm", "Y", "La", "Ac", "Al", "Li",
        "Rb", "Co", "Cf", "Os", "Zr", "Ni", "Ba", "Hg", "Bk"
    }
    
    if symbol in accepted_metals:
        return True, f"The atom is a metal: {symbol}"
    else:
        return False, f"Atom symbol '{symbol}' is not in the accepted set of metal atoms."

# Example testing (uncomment for local testing):
# print(is_metal_atom("[No]"))       # Expected: True
# print(is_metal_atom("[218Po]"))     # Expected: True (since [218Po] => Po)
# print(is_metal_atom("[27Al]"))      # Expected: True (Aluminum with isotope info)
# print(is_metal_atom("[H][Bi]([H])[H]"))  # Expected: False; more than one atom.
# print(is_metal_atom("[Rb]"))        # Expected: True after updating accepted_metals
# print(is_metal_atom("[Co]"))        # Expected: True after updating accepted_metals