"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: Mineral (inorganic, often ionic/geologically formed species)
Based on revised heuristics:
  1. Minerals are often ionic salts – their SMILES will usually consist of multiple disconnected fragments.
  2. Even if a molecule is a single fragment, the presence of metals/non‐organic elements (elements not common in organic compounds)
     is a strong indicator of mineral identity.
  3. Organic molecules (even if they carry formal charges, as in fatty acids) are generally a single covalently bound entity with only
     the “organic set” elements: H, C, N, O, F, P, S, Cl, Br, I.
Examples include calcium difluoride, greigite, various hydrates and ionic salts.
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is considered a mineral based on its SMILES string.
    Revised heuristic criteria:
      1. Parse the SMILES string.
      2. Count disconnected fragments; multiple fragments are typical for ionic mineral salts.
      3. If there is only a single fragment, check for the presence of any metal or non‐organic element.
         We define the organic set as: H, C, N, O, F, P, S, Cl, Br, and I.
         The presence of an element outside this set (e.g., Na, Ca, Fe, etc.) strongly favors a mineral.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets our criteria for a mineral, False otherwise.
        str: A reason describing the classification.
    """
    # Try to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get disconnected fragments in the molecule.
    frags = Chem.GetMolFrags(mol, asMols=True)
    frag_count = len(frags)
    
    # Define the atomic numbers of elements common in organic compounds:
    # H (1), C (6), N (7), O (8), F (9), P (15), S (16), Cl (17), Br (35), I (53)
    organic_set = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}
    
    # Check if any atom in the entire molecule is not in the organic set.
    metal_present = any(atom.GetAtomicNum() not in organic_set for atom in mol.GetAtoms())
    
    # Primary criteria:
    # 1. If there are multiple fragments, that is typical for ionic, geologically formed salts.
    if frag_count > 1:
        return True, f"Multiple fragments detected (count: {frag_count}), typical for ionic mineral compounds."
    
    # 2. For single-fragment species, if there is any element that is not typical for organic compounds,
    #    then it is likely a mineral (or an inorganic species).
    if frag_count == 1 and metal_present:
        return True, "Single fragment but contains element(s) not typical in organic compounds, indicating possible mineral nature."
    
    # Otherwise assume it is not a mineral.
    return False, "Does not display inorganic or ionic features typical of minerals."

# Example calls for testing (uncomment to test):
# print(is_mineral("[S--].[S--].[S--].[S--].[Fe++].[Fe+3].[Fe+3]"))  # greigite
# print(is_mineral("[F-].[F-].[Ca++]"))  # calcium difluoride
# print(is_mineral("CCCCCCCCCCCCCCCCCC(O)C([O-])=O"))  # 2-hydroxyarachidate (organic acid) should not be classified as mineral