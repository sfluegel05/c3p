"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: A derivative of the dimethylisoalloxazine (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton with a substituent on the 10 position (i.e. flavin).
  
This program uses a SMARTS pattern that approximates the isoalloxazine core. The pattern
“c1cc2nc3c(nc(=O)[nH]c3=O)n([*:1])c2cc1” is designed to:
  - Match a fused, three‐ring heterocycle corresponding approximately to the isoalloxazine skeleton.
  - Require (via n([*:1])) that one of the ring nitrogens (position 10) bears a substituent (which can be any atom or group).
  
This is only one possible approach to classify flavins.
"""

from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin derivative (a derivative of the dimethylisoalloxazine skeleton with a substituent at the 10 position)
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a flavin derivative, False otherwise.
        str: Reason for the classification decision.
    """
    # Attempt to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the isoalloxazine core with a mandatory substituent on the bridging N.
    # The SMARTS “c1cc2nc3c(nc(=O)[nH]c3=O)n([*:1])c2cc1” is based on the skeleton of dimethylisoalloxazine.
    # It requires that one of the heterocyclic nitrogens carries a substituent (the [*:1] part) which is assumed to be at the 10 position.
    isoalloxazine_smarts = "c1cc2nc3c(nc(=O)[nH]c3=O)n([*:1])c2cc1"
    core = Chem.MolFromSmarts(isoalloxazine_smarts)
    if core is None:
        # In case of an error building the pattern, we cannot proceed.
        return False, "Error constructing isoalloxazine SMARTS pattern"
    
    # Check if the molecule contains the isoalloxazine core with a substituent on the bridging N.
    if not mol.HasSubstructMatch(core):
        return False, "Does not contain the required isoalloxazine core with a substituent at the 10 position"
    
    return True, "Contains isoalloxazine core (dimethylisoalloxazine skeleton) with a substituent at N10, classified as flavin"

# Example test cases (uncomment to run)
# test_smiles = [
#     "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C",  # lumiflavin example
#     "CC1=C(C)C=C2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3=NC2=C1",  # riboflavin example
#     "CCCC",  # non-flavin
# ]
# for smi in test_smiles:
#     result, reason = is_flavin(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")