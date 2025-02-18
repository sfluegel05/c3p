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
    
    The function first tries to match a contiguous substructure motif defined as:
      [NX3;!$(N-C(=O))]-[C;H1,H2]-[C](=O)-[O]-[CX4]
    That is, a free (unacylated) amino group attached to an α–carbon (which is CH or CH2)
    attached to a carbonyl that is esterified (–O– connected to an sp³ carbon).
    
    Post–processing is done to ensure the match is chemically plausible:
      • The α–carbon (the second atom matched) must have at least 1 hydrogen.
      • The oxygen (fourth atom in the fragment) must be bonded to a carbon (other than
        the carbonyl carbon) that is sp³‐hybridized.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one alpha–amino acid ester motif is found, False otherwise.
        str: Explanation for the classification.
    """
    # Try to parse the input SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that requires the contiguous fragment:
    #     [NX3;!$(N-C(=O))]   -> free (non acylated) amino nitrogen
    #     -[C;H1,H2]         -> an α–carbon carrying one or two hydrogens
    #     -[C](=O)           -> a carbonyl carbon
    #     -[O]               -> the carbonyl oxygen (in an ester)
    #     -[CX4]             -> an sp³ carbon (the ester alkyl group)
    aa_ester_smarts = "[NX3;!$(N-C(=O))]-[C;H1,H2]-[C](=O)-[O]-[CX4]"
    aa_ester_pattern = Chem.MolFromSmarts(aa_ester_smarts)
    if aa_ester_pattern is None:
        return False, "Error in SMARTS pattern construction"
    
    # Look for all substructure matches of the motif.
    matches = mol.GetSubstructMatches(aa_ester_pattern)
    if not matches:
        return False, "No alpha–amino acid ester motif found"
    
    # For each match, perform additional checks.
    # The match tuple corresponds to atoms (in order): [amino N, α–C, carbonyl C, ester O, alkyl C]
    valid_match_found = False
    for match in matches:
        amino_idx, alpha_idx, carbonyl_idx, esterO_idx, ester_alkyl_idx = match
        
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        # Check that the alpha carbon has at least one hydrogen 
        # (the SMARTS [C;H1,H2] should enforce this, but we double-check).
        if alpha_atom.GetTotalNumHs() < 1:
            continue

        # Check that the ester oxygen (atom at esterO_idx) has a neighbor (besides the carbonyl)
        # that is an sp³ carbon. (We already demanded it via SMARTS; but in case of conjugated systems,
        # this extra check helps to avoid mis‐assignment in larger scaffolds.)
        esterO_atom = mol.GetAtomWithIdx(esterO_idx)
        valid_neighbor = False
        for nbr in esterO_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_idx:
                continue  # skip the carbonyl carbon (which is already in the match)
            # If the neighbor is an sp3 carbon, we count it as valid.
            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.HybridizationType.SP3:
                valid_neighbor = True
                break
        if not valid_neighbor:
            continue

        # If we pass all additional checks, then we confidently classify as an alpha–amino acid ester.
        valid_match_found = True
        break

    if valid_match_found:
        return True, "Contains an alpha–amino acid ester moiety"
    else:
        return False, "No valid alpha–amino acid ester motif found after additional checks"

# Example use (for testing purposes - uncomment as needed):
# test_smiles = [
#     "O(C(=O)C(N)C(OCC)=O)CC",           # Diethyl aminomalonate - should be True
#     "OC(=O)CN",                        # methyl glycinate - should be True
#     "CCCC",                            # Negative test - should be False
# ]
# for smi in test_smiles:
#     result, reason = is_alpha_amino_acid_ester(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")