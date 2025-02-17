"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: A metal atom.
A metal atom is defined as an atom that exhibits typical metallic properties,
being typically shiny with high electrical and thermal conductivity.
For this classifier, the SMILES string must represent exactly a single neutral atom
(with any isotopic labeling removed). The atomic symbol must be in the accepted set.
The accepted set of metal atoms has been expanded following feedback.
"""

from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    
    Criteria:
      - The SMILES must parse to a molecule.
      - The molecule must be composed of a single disconnected fragment.
      - That fragment must contain exactly one atom.
      - The atom must have a formal charge of 0.
      - The atomic symbol (ignoring any isotope numbers) must be in the accepted set of metal atoms.
    
    Accepted metal symbols include (expanded): 
      No, Cn, Po, Cr, Pa, Tc, Cs, Sr, Bi, Tb, Bh, Sn, Gd, Eu, Fm, Mo, Sc, Ho, Md, Mn,
      K, Am, Sm, Ce, Nb, Tm, Y, La, Ac, Al, Li, Rb, Co, Cf, Os, Zr, Ni, Ba, Hg, Bk,
      Mt, Ru, Dy, Cd, Rh, Zn, Cm, Er, Ra, Ti.
      
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
    
    # Get disconnected fragments from the molecule.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) != 1:
        return False, f"Expected a single fragment, but got {len(frags)} fragments."
    
    frag = frags[0]
    # Ensure the fragment contains exactly one atom.
    if frag.GetNumAtoms() != 1:
        return False, f"Expected 1 atom, but the fragment contains {frag.GetNumAtoms()} atoms."
    
    # Get the single atom and ensure it is neutral (formal charge 0).
    atom = frag.GetAtomWithIdx(0)
    if atom.GetFormalCharge() != 0:
        return False, f"The atom has a formal charge of {atom.GetFormalCharge()}, but a neutral atom is expected."
    
    # Obtain the atomic symbol (this ignores isotope information)
    symbol = atom.GetSymbol()
    
    # Define the accepted set of metal atoms (expanded based on feedback)
    accepted_metals = {
        "No", "Cn", "Po", "Cr", "Pa", "Tc", "Cs", "Sr", "Bi", "Tb", "Bh",
        "Sn", "Gd", "Eu", "Fm", "Mo", "Sc", "Ho", "Md", "Mn", "K",
        "Am", "Sm", "Ce", "Nb", "Tm", "Y", "La", "Ac", "Al", "Li",
        "Rb", "Co", "Cf", "Os", "Zr", "Ni", "Ba", "Hg", "Bk",
        "Mt", "Ru", "Dy", "Cd", "Rh", "Zn", "Cm", "Er", "Ra", "Ti"
    }
    
    if symbol in accepted_metals:
        return True, f"The atom is a metal: {symbol}"
    else:
        return False, f"Atom symbol '{symbol}' is not in the accepted set of metal atoms."

# Example testing (uncomment for local testing):
# print(is_metal_atom("[No]"))       # Expected: True
# print(is_metal_atom("[Cn]"))       # Expected: True
# print(is_metal_atom("[218Po]"))     # Expected: True (interprets as 'Po')
# print(is_metal_atom("[27Al]"))      # Expected: True (Aluminum from isotope annotation)
# print(is_metal_atom("[H][Bi]([H])[H]"))  # Expected: False; more than one atom.
# print(is_metal_atom("[Rb]"))        # Expected: True
# print(is_metal_atom("[Co]"))        # Expected: True