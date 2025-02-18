"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: alpha-amino acid ester
Definition: The amino acid ester derivative obtained by the formal condensation 
of an alpha–amino acid with an alcohol.
An α–amino acid typically has the structure H2N–CHR–COOH. Upon esterification, 
the –COOH is converted to –COOR. Here we search for a contiguous motif corresponding to:
   N–αC–C(=O)–O–C
with the refinements that (a) the amino nitrogen is not acylated,
(b) the α–carbon is either CH or CH2 (ensuring at least one hydrogen),
and (c) the ester oxygen is bound to an sp³ carbon.
"""

from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha–amino acid ester based on its SMILES string.
    
    The function looks for a contiguous fragment defined as:
      [NX3;!$(N-C(=O))]-[C;H1,H2]-[C](=O)-[O]-[CX4]
    That is, a free (non-acylated) amino nitrogen attached to an α–carbon (which is CH or CH2)
    attached to a carbonyl that is esterified (–O– attached to an sp³ carbon).
    
    Additional checks are performed to verify:
      • the α–carbon has at least one hydrogen,
      • the ester oxygen (in the carbonyl ester group) is bound to an sp³ carbon (other than the carbonyl).
    
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
    
    # Define the SMARTS pattern for the alpha-amino acid ester fragment.
    aa_ester_smarts = "[NX3;!$(N-C(=O))]-[C;H1,H2]-[C](=O)-[O]-[CX4]"
    aa_ester_pattern = Chem.MolFromSmarts(aa_ester_smarts)
    if aa_ester_pattern is None:
        return False, "Error constructing SMARTS pattern"
    
    # Find all matches of the pattern.
    matches = mol.GetSubstructMatches(aa_ester_pattern)
    if not matches:
        return False, "No alpha–amino acid ester motif found"
    
    # Check each match for additional chemical plausibility.
    valid_match_found = False
    for match in matches:
        # In some cases the match tuple may have extra indices; we only require the first 5.
        if len(match) < 5:
            continue
        amino_idx, alpha_idx, carbonyl_idx, esterO_idx, ester_alkyl_idx = match[:5]
        
        # Check that the α–carbon has at least one hydrogen.
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        if alpha_atom.GetTotalNumHs() < 1:
            continue
        
        # Verify that the ester oxygen (esterO_idx) is connected to an sp³ carbon other than the carbonyl.
        esterO_atom = mol.GetAtomWithIdx(esterO_idx)
        valid_neighbor = False
        for nbr in esterO_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_idx:
                continue  # Skip the carbonyl carbon.
            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.HybridizationType.SP3:
                valid_neighbor = True
                break
        if not valid_neighbor:
            continue
        
        # If all check are passed, we have a valid match.
        valid_match_found = True
        break

    if valid_match_found:
        return True, "Contains a valid alpha–amino acid ester moiety"
    else:
        return False, "No valid alpha–amino acid ester motif found after additional checks"

# Example use (for testing):
if __name__ == "__main__":
    test_smiles = [
        "O(C(=O)C(N)C(OCC)=O)",  # Diethyl aminomalonate, expected True
        "OC(=O)CN",             # methyl glycinate, expected True
        "CCCC",                 # Negative test, expected False
    ]
    for smi in test_smiles:
        result, reason = is_alpha_amino_acid_ester(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")