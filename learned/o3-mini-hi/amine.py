"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: Amine
Definition: A compound formally derived from ammonia by replacing one, two or three hydrogen atoms
by hydrocarbyl groups. Valid amines include primary, secondary and tertiary amines. (A nitrogen directly
bonded to a carbonyl - as in an amide - is excluded.)
"""

from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule contains an amine functional group based on a set of SMARTS patterns.
    
    A valid amine is (i) primary [–NH2], (ii) secondary [–NH–] or (iii) tertiary [–N–] provided that 
    the nitrogen is not directly attached to a carbonyl carbon (i.e. not part of an amide) and is not quaternary.
    
    The approach uses SMARTS patterns for:
      - primary amine: "[NX3;H2;!$(N-C=O)]"
      - secondary amine: "[NX3;H1;!$(N-C=O)]"
      - tertiary amine: "[NX3;H0;!$(N-C=O);!$([N+])]"
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one qualifying amine group is found, False otherwise.
              (Note: For very complex molecules the result might be ambiguous.)
        str: An explanation or reason for the classification.
    """
    
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for primary, secondary and tertiary amines.
    # The substructure [NX3] selects trivalent N; Hx selects the number of attached hydrogens.
    # The negative SMARTS !$(N-C=O) prevents matching when the nitrogen is directly
    # bonded to a carbonyl carbon (as in an amide). For tertiary amines we also exclude positive charges.
    amine_smarts = {
        "primary": "[NX3;H2;!$(N-C=O)]",
        "secondary": "[NX3;H1;!$(N-C=O)]",
        "tertiary": "[NX3;H0;!$(N-C=O);!$([N+])]"
    }
    
    # Check each pattern in turn
    for typ, smarts in amine_smarts.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue  # skip pattern if invalid (should not happen)
        if mol.HasSubstructMatch(patt):
            return True, f"Found a {typ} amine substructure that is not bound to a carbonyl."
    
    # If no pattern matched, no qualifying amine was detected.
    return False, "No qualifying amine functional group found."

# Example usage:
# result, reason = is_amine("CNc1ccccc1")
# print(result, reason)