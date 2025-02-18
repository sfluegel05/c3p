"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: alpha–amino acid ester
Definition: The amino acid ester derivative obtained by the formal condensation 
of an alpha–amino acid with an alcohol.
A typical alpha–amino acid is H2N–CHR–COOH. In an ester derivative the –COOH is converted 
to –COOR. 

This program looks for an uninterrupted fragment:
    N – alpha C – C(=O) – O – Alkyl
with the following refinements:
  (a) The amino nitrogen must not be acylated.
  (b) The α–carbon must possess at least one hydrogen (checked after matching).
  (c) The ester oxygen must be connected to an sp³ carbon other than the carbonyl.

If a match is found that passes the additional checks, then the structure is classified 
as an alpha–amino acid ester.
"""

from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha–amino acid ester based on its SMILES string.
    
    This function searches for the contiguous fragment:
      [NX3;!$(N-C(=O))]-[C]-C(=O)-O-[CX4]
    which represents a free (non-acylated) amino nitrogen attached to an α–carbon,
    attached to a carbonyl that is esterified.
    
    Additional checks are applied:
      • the α–carbon has at least one hydrogen (using GetTotalNumHs), and 
      • the ester oxygen is connected to an sp³ carbon (other than the carbonyl carbon).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a valid alpha–amino acid ester motif, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that is less strict about explicit hydrogen counts on the α–carbon.
    # The pattern:
    #   [NX3;!$(N-C(=O))]  --> free amino nitrogen.
    #   -[C]               --> any carbon (we will later ensure it has at least one hydrogen).
    #   -C(=O)             --> carbonyl carbon.
    #   -O-                --> ester oxygen.
    #   [CX4]              --> an sp3 carbon (alkyl group).
    aa_ester_smarts = "[NX3;!$(N-C(=O))]-[C]-C(=O)-O-[CX4]"
    aa_ester_pattern = Chem.MolFromSmarts(aa_ester_smarts)
    if aa_ester_pattern is None:
        return False, "Error constructing SMARTS pattern"
    
    # Find all substructure matches.
    matches = mol.GetSubstructMatches(aa_ester_pattern)
    if not matches:
        return False, "No alpha–amino acid ester motif found"
    
    # Check additional properties for each found match.
    for match in matches:
        # Expecting a 5-atom match: (amino N, α–carbon, carbonyl C, ester O, alkyl C)
        if len(match) < 5:
            continue  # Skip incomplete matches
        
        amino_idx, alpha_idx, carbonyl_idx, esterO_idx, ester_alkyl_idx = match[:5]
        
        # Check that the α–carbon has at least one hydrogen. Implicit hydrogens count.
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        if alpha_atom.GetTotalNumHs() < 1:
            continue
        
        # Verify that the ester oxygen is connected to an sp³ carbon that is not the carbonyl.
        esterO_atom = mol.GetAtomWithIdx(esterO_idx)
        valid_neighbor = False
        for nbr in esterO_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_idx:
                continue  # Skip the carbonyl atom.
            # Check that the neighbor is carbon and is sp3 hybridized.
            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.HybridizationType.SP3:
                valid_neighbor = True
                break
        if not valid_neighbor:
            continue
        
        # Found a matching fragment that passes the checks.
        return True, "Contains a valid alpha–amino acid ester moiety"
    
    # If no match passed the additional checks, return False.
    return False, "No valid alpha–amino acid ester motif found after additional checks"

# Example testing (uncomment the block below if running as main)
if __name__ == "__main__":
    test_examples = [
        # Expected True examples:
        "O(C(=O)C(N)C(OCC)=O)",  # Diethyl aminomalonate
        "O(CC1=CC=CC=C1)C(=O)CN", # Benzyl glycinate
        "OC(=O)CN",              # methyl glycinate
        # Expected False example:
        "CCCC"
    ]
    for smi in test_examples:
        result, reason = is_alpha_amino_acid_ester(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")