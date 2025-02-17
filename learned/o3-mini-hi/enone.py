"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: enone – alpha,beta-unsaturated ketones
Definition: An enone is an alpha,beta-unsaturated ketone of general formula 
            R(1)R(2)C=CR(3)-C(=O)R(4) (with R(4) ≠ H) in which the carbonyl group (C=O)
            is directly conjugated to a C=C double bond.
            
Algorithm:
  1. Parse the input SMILES string.
  2. Search for a fragment having the pattern [#6]=[#6]-[C](=O).
     This pattern identifies any carbon (atom1) double-bonded to a second carbon (atom2)
     that is single-bonded to a carbonyl group (atom3, with =O).
  3. Post-filter each hit by checking that the carbonyl carbon (atom3) does not have any hydrogens attached.
     (This step prevents aldehydes from being classified as enones.)
  4. Return True with a reason if a valid match is found; otherwise False.
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone (an alpha,beta-unsaturated ketone)
    based on its SMILES string.
    
    The algorithm first searches for a C=C–C(=O) fragment via a SMARTS query.
    Then, it ensures that the carbonyl carbon carries no hydrogen atoms, which excludes aldehydes.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a valid enone motif is found, otherwise False.
        str: A message explaining the classification result.
    """
    
    # Parse the SMILES string to a RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that searches for a potential enone fragment:
    #   [#6]=[#6]-[C](=O)
    # Explanation:
    #   [#6]         : any carbon atom,
    #   =[#6]        : double-bonded to another carbon,
    #   -[C](=O)    : single bond to a carbon that is double-bonded to an oxygen.
    pattern = Chem.MolFromSmarts("[#6]=[#6]-[C](=O)")
    if pattern is None:
        return None, None  # fallback if SMARTS creation fails (should not happen)
    
    matches = mol.GetSubstructMatches(pattern)
    # Evaluate each match to see if it qualifies as an enone.
    for match in matches:
        # match is a tuple: (index_of_first_C, index_of_second_C, index_of_carbonyl_C)
        carbonyl_atom = mol.GetAtomWithIdx(match[2])
        # Check that the carbonyl carbon has no attached hydrogens (i.e. not an aldehyde)
        if carbonyl_atom.GetTotalNumHs() == 0:
            reason = ("Contains a conjugated enone motif (C=C–C(=O)) found at atom indices: "
                      "C1 (idx {}), C2 (idx {}), and carbonyl carbon (idx {}).".format(match[0], match[1], match[2]))
            return True, reason
    
    # No valid enone fragment found so far.
    return False, "Does not contain a conjugated enone motif (no valid C=C–C(=O) fragment with non-aldehyde carbonyl found)"

# For basic testing you can uncomment the lines below:
# test_smiles = [
#     "C=CC(=O)C",                # minimal enone (methyl vinyl ketone derivative)
#     "CC1=CC(=O)CC(C)(C)C1",      # isophorone, which qualifies as an enone
#     "O=C(C)C=C",                # an enal (aldehyde) – should not be classified as enone
# ]
# for sm in test_smiles:
#     result, msg = is_enone(sm)
#     print("SMILES:", sm)
#     print("Result:", result)
#     print("Reason:", msg, "\n")