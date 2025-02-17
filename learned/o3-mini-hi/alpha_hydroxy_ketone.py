"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: Alpha-hydroxy ketone
Definition: A ketone containing a hydroxy group on the alpha-carbon relative to the C=O group.
Improvement notes:
 - Instead of iterating over ketone matches and then checking neighbors,
   we combine the motif into a SMARTS pattern.
 - We require a C=O flanked by two carbons, with one of those carbons bearing an –OH.
 - Two SMARTS patterns are used (one for an alpha–OH on the left and one for the right).
 - We add explicit hydrogens so that –OH groups are recognized.
 - (Further filtering might be needed to exclude cases such as highly decorated sugars.)
"""

from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined as a ketone (C=O with two carbon substituents)
    that has at least one hydroxy (OH) group directly bonded to one of the carbons adjacent
    to the carbonyl group.
    
    Two complementary SMARTS patterns are used:
      Pattern1: matches structures of the form: R–C(OH)–C(=O)–R'
      Pattern2: matches structures of the form: R–C(=O)–C(OH)–R
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: A reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydroxyl groups appear explicitly.
    mol = Chem.AddHs(mol)
    
    # Define SMARTS patterns.
    # Pattern1 looks for an aliphatic carbon (any, marked as [#6]) attached to 
    # a carbon that bears a directly attached –OH (represented as -[OX2H]),
    # followed by a carbonyl group (C(=O)) which in turn is attached to another carbon.
    # It effectively encodes:  R – C(-OH) – C(=O) – R'
    pattern1 = Chem.MolFromSmarts("[#6]!@[$([C](-[OX2H]))]!@C(=O)[#6]")
    
    # Pattern2 is essentially the mirror image:  R – C(=O) – C(-OH) – R
    pattern2 = Chem.MolFromSmarts("[#6]C(=O)[!@][$([C](-[OX2H]))]!@[#6]")
    
    # Check if either pattern is found in the molecule.
    if mol.HasSubstructMatch(pattern1):
        return True, "Found a ketone group with an alpha hydroxy substituent (pattern1)"
    if mol.HasSubstructMatch(pattern2):
        return True, "Found a ketone group with an alpha hydroxy substituent (pattern2)"
    
    # If no match is found, return false.
    return False, "No substructure matching alpha-hydroxy ketone pattern found"

# Example usage:
if __name__ == "__main__":
    # This is one of the known positives: (3S)-1-hydroxy-3-methylpentan-2-one
    test_smiles = "OCC([C@H](CC)C)=O"
    result, reason = is_alpha_hydroxy_ketone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)