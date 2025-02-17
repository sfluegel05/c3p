"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: alpha-amino acid ester
Definition: The amino acid ester derivative obtained by the formal condensation 
of an alpha–amino acid with an alcohol.
An α–amino acid typically has the structure H2N–CHR–COOH.
In an ester derivative, the –COOH is converted to –COOR.
We search for a contiguous fragment corresponding to:
   N–αC–C(=O)–O–C
with the refinements that (a) the amino nitrogen is not acylated,
(b) the α–carbon must have at least one hydrogen (i.e. be CH or CH2, including chiral variants), and 
(c) the ester oxygen is bound to an sp³ carbon (other than the carbonyl carbon).
"""

from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha–amino acid ester based on its SMILES string.
    
    This function looks for a contiguous fragment defined as:
      [NX3;!$(N-C(=O))]-[CH1,CH2]-C(=O)-O-[CX4]
    i.e. a free (non-acylated) amino nitrogen attached to an α–carbon 
    (which is CH or CH2, including chiral forms) attached to a carbonyl 
    that is esterified (-O- attached to an sp³ carbon).
    
    Additional checks:
      • the α–carbon has at least one hydrogen (ensured by SMARTS but double‐checked),
      • the ester oxygen is connected to an sp³ carbon (other than the attached carbonyl).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a valid alpha–amino acid ester motif is found, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more flexible SMARTS pattern.
    # The pattern now is:
    #   [NX3;!$(N-C(=O))] : a free (non-acylated) amino nitrogen
    #   -[CH1,CH2]       : an alpha-carbon with 1 or 2 hydrogens (chiral centers included)
    #   -C(=O)           : a carbonyl carbon
    #   -O-              : an ester oxygen
    #   [CX4]            : bound to an sp3 carbon (alkyl group)
    aa_ester_smarts = "[NX3;!$(N-C(=O))]-[CH1,CH2]-C(=O)-O-[CX4]"
    aa_ester_pattern = Chem.MolFromSmarts(aa_ester_smarts)
    if aa_ester_pattern is None:
        return False, "Error constructing SMARTS pattern"
    
    # Find all matches of the pattern.
    matches = mol.GetSubstructMatches(aa_ester_pattern)
    if not matches:
        return False, "No alpha–amino acid ester motif found"
    
    # Check additional properties in each match.
    for match in matches:
        # Expect match to be a 5-atom tuple: (amino N, alpha carbon, carbonyl C, ester O, alkyl C)
        if len(match) < 5:
            continue  # Skip incomplete matches
        
        amino_idx, alpha_idx, carbonyl_idx, esterO_idx, ester_alkyl_idx = match[:5]
        
        # Check that the alpha-carbon has at least one hydrogen.
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        if alpha_atom.GetTotalNumHs() < 1:
            continue
        
        # Verify that the ester oxygen (esterO_idx) is connected to an sp³ carbon 
        # that is not the carbonyl carbon.
        esterO_atom = mol.GetAtomWithIdx(esterO_idx)
        valid_neighbor = False
        for nbr in esterO_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_idx:
                continue  # Skip the carbonyl carbon.
            # Check that neighbor is carbon and is sp3.
            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.HybridizationType.SP3:
                valid_neighbor = True
                break
        if not valid_neighbor:
            continue
        
        # If all conditions are satisfied, return True.
        return True, "Contains a valid alpha–amino acid ester moiety"
    
    # If no match passes additional checks, then return False.
    return False, "No valid alpha–amino acid ester motif found after additional checks"

# For testing (example use):
if __name__ == "__main__":
    test_smiles = [
        "O(C(=O)C(N)C(OCC)=O)",  # Diethyl aminomalonate (expected True)
        "O(CC1=CC=CC=C1)C(=O)CN", # Benzyl glycinate (expected True)
        "OC(=O)CN",              # methyl glycinate (expected True)
        "CCCC",                  # Negative test, expected False
    ]
    for smi in test_smiles:
        result, reason = is_alpha_amino_acid_ester(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")