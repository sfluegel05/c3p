"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: enone – alpha,beta-unsaturated ketones
Definition: A molecule is an enone if it contains a C=C–C(=O) fragment in which
            the carbonyl carbon is not an aldehyde (i.e. it carries no hydrogen).
This SMARTS-based approach looks for a carbon–carbon double bond directly connected 
via a single bond to a carbonyl carbon (which must be substituted by a carbon).
The used SMARTS is: [#6]=[#6]-[C;!H1](=O)
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone (an alpha,beta-unsaturated ketone)
    based on its SMILES string.
    
    The algorithm uses a SMARTS pattern to search for the defining fragment:
      – A carbon–carbon double bond,
      – Followed by a single bond to a carbonyl carbon (C(=O))
      – In which the carbonyl carbon does not have any attached hydrogen (not an aldehyde).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if an enone motif is found, otherwise False.
        str: A reason explaining the classification.
             If no valid motif is found the reason explains the missing fragment.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for the enone fragment:
    # Explanation:
    #   [#6]         : any carbon atom.
    #   =[#6]        : double-bonded to a carbon atom.
    #   -[C;!H1]    : connected via a single bond to a carbon that does NOT have any attached hydrogen
    #                  (i.e. to avoid aldehydes).
    #   (=O)        : that carbon is double-bonded to an oxygen.
    enone_smarts = "[#6]=[#6]-[C;!H1](=O)"
    enone_pat = Chem.MolFromSmarts(enone_smarts)
    if enone_pat is None:
        return None, None  # if SMARTS pattern somehow fails

    # Search for the enone substructure in the molecule.
    matches = mol.GetSubstructMatches(enone_pat)
    if matches:
        # Get first match and report the indices of the atoms involved.
        # The match is a tuple of three indices corresponding to:
        #   idx0: carbon involved in the C=C bond (first carbon)
        #   idx1: second carbon of the C=C bond (alpha carbon)
        #   idx2: carbonyl carbon
        match = matches[0]
        reason = ("Contains a conjugated enone motif (C=C–C(=O)) found at atom indices: "
                  "C1 (idx {}), C2 (idx {}), and carbonyl carbon (idx {}).".format(match[0], match[1], match[2]))
        return True, reason
    else:
        return False, "Does not contain a conjugated enone motif (no C=C–C(=O) fragment was found)"

# For basic testing you can uncomment the lines below:
# test_smiles = [
#     "C=CC(=O)C",                # minimal enone (methyl vinyl ketone derivative)
#     "CC1=CC(=O)CC(C)(C)C1",      # isophorone
#     "O=C(C)C=C",                # an enal (aldehyde) – should not be classified as enone
# ]
# for sm in test_smiles:
#     result, msg = is_enone(sm)
#     print("SMILES:", sm)
#     print("Result:", result)
#     print("Reason:", msg, "\n")