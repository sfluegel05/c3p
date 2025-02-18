"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: beta-lactam
A beta-lactam is defined as a lactam in which the amide bond is contained
within a four-membered ring that includes the amide nitrogen and the carbonyl carbon.
We use a SMARTS pattern that matches a 4-membered ring with a nitrogen and a carbonyl.
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    
    We first sanitize the input molecule. Then, to avoid misclassification due to salts
    or multiple fragments, we select the largest fragment. Finally, we use a SMARTS pattern
    that matches a 4-membered ring (using the "R4" qualifier) which contains a nitrogen
    atom and a carbon atom that is double-bonded to an oxygen.
    
    SMARTS explanation:
      - [NX3;R4]     => any trivalent nitrogen that is a member of a ring of exactly 4 atoms.
      - [C;R4](=O)   => any carbon in a ring of 4 that is double-bonded to an oxygen.
      - [C;R4]      => a 4-membered ring carbon (the remaining ring atom).
      - [C;R4]      => the fourth atom in the ring.
    The ring order in the SMARTS pattern can match in any rotational order.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a beta-lactam, False otherwise.
        str: A reason supporting the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove salts by selecting the largest fragment
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments could be parsed from SMILES"
    # Pick the fragment with the most heavy atoms as the main fragment
    main_frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # Define a SMARTS pattern for a beta-lactam ring.
    # This pattern looks for a 4-membered ring ("R4") that contains:
    # - A nitrogen atom ([NX3;R4])
    # - A carbon that is double-bonded to an oxygen ([C;R4](=O))
    # - And two additional ring atoms ([C;R4] twice)
    beta_lactam_smarts = "[NX3;R4][C;R4](=O)[C;R4][C;R4]"
    pattern = Chem.MolFromSmarts(beta_lactam_smarts)
    if pattern is None:
        return False, "Invalid beta-lactam SMARTS pattern"
    
    # Search for the beta-lactam motif in the main fragment
    if main_frag.HasSubstructMatch(pattern):
        return True, "Beta-lactam ring detected: 4-membered ring with a carbonyl carbon and an adjacent nitrogen"
    
    return False, "No beta-lactam ring found (4-membered ring with the required amide bond was not detected)"

# Example usage:
if __name__ == "__main__":
    # Testing on the simplest beta-lactam: azetidin-2-one
    test_smiles = "O=C1CCN1"
    result, reason = is_beta_lactam(test_smiles)
    print("Result:", result)
    print("Reason:", reason)