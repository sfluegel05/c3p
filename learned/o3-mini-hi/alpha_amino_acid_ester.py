"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: alpha-amino acid ester
Definition: The amino acid ester derivative obtained by the formal condensation 
of an alpha-amino acid with an alcohol.
An α‐amino acid typically has the structure H2N–CHR–COOH. Upon esterification, 
the –COOH is converted to –COOR. Here we search for a substructure corresponding to:
   N–αC–C(=O)O[*] 
with the refinements that (a) the nitrogen is not acylated (i.e. not an amide) and 
(b) the α‐carbon is either CH (one hydrogen) or CH2 (two hydrogens) so that glycine and 
substituted amino acids are both covered.
"""

from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    This function looks for the motif:
        [NX3;!$(N-C(=O))]-[C;H1,H2]-[C](=O)O[*]
    which represents a (free) amino group attached to an alpha carbon (with one or two H's)
    that is further bound to a carbonyl carbon in an ester (C(=O)O-).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains an alpha-amino acid ester motif, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an improved SMARTS pattern:
    # - [NX3;!$(N-C(=O))] ensures the amino nitrogen is not immediately acylated.
    # - [C;H1,H2] requires the alpha carbon be either CH or CH2.
    # - [C](=O)O[*] looks for a carbonyl group that is esterified.
    aa_ester_smarts = "[NX3;!$(N-C(=O))]-[C;H1,H2]-[C](=O)O[*]"
    aa_ester_pattern = Chem.MolFromSmarts(aa_ester_smarts)
    if aa_ester_pattern is None:
        return False, "Error in SMARTS pattern"

    # Check if the molecule has at least one match to the defined alpha-amino acid ester motif.
    if not mol.HasSubstructMatch(aa_ester_pattern):
        return False, "No alpha-amino acid ester motif found"

    return True, "Contains an alpha-amino acid ester moiety"

# Example use (uncomment for testing):
# smiles_examples = [
#     "O(C(=O)C(N)C(OCC)=O)CC",      # Diethyl aminomalonate
#     "OC(=O)CN",                   # methyl glycinate
#     "CCCC",                       # negative test
# ]
# for smi in smiles_examples:
#     result, reason = is_alpha_amino_acid_ester(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")