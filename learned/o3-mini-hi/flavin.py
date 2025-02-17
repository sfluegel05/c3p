"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: A derivative of the dimethylisoalloxazine 
(7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton 
with a substituent on the 10 position (i.e. flavin).

This program uses a SMARTS pattern that requires:
  - An aromatic benzene ring with two explicit methyl groups 
    (at positions corresponding to 7 and 8), and
  - A fused heterocyclic system (the isoalloxazine core) where 
    a bridging nitrogen carries a non-hydrogen substituent.
Only molecules matching this assembly are classified as flavin.
"""

from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin derivative based on its SMILES string.
    A flavin is defined as a derivative of the dimethylisoalloxazine (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione)
    skeleton with a substituent on the 10 position. In our approach we enforce that the skeleton contains 
    two methyl groups on the benzene ring and that the bridging nitrogen (N10) has a non-hydrogen substituent.

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        bool: True if the molecule is classified as a flavin derivative, False otherwise.
        str: A string explaining the reason for the classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the dimethylisoalloxazine core with the N10 substituent.
    # The pattern "c1c(C)cc2nc3c(nc(=O)[nH]c3=O)n([*:1])c2c1C" is designed to:
    #  - Begin with an aromatic ring (c1...) where the first atom has an attached methyl (c(C)).
    #  - Fuse into a heterocyclic system representing the isoalloxazine core.
    #  - Require that one of the nitrogens (the bridging, N10) bears a substituent (n([*:1])).
    #  - End the benzene ring with an explicit methyl (c1C) to enforce the dimethyl pattern.
    isoalloxazine_smarts = "c1c(C)cc2nc3c(nc(=O)[nH]c3=O)n([*:1])c2c1C"
    core = Chem.MolFromSmarts(isoalloxazine_smarts)
    
    if core is None:
        return False, "Error constructing isoalloxazine SMARTS pattern"
    
    # Check if the molecule contains the isoalloxazine core according to our definition.
    if not mol.HasSubstructMatch(core):
        return False, "Does not contain the required dimethylisoalloxazine core with a substituent at N10"
    
    return True, "Contains isoalloxazine core (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) with a substituent at N10, classified as flavin"

# Uncomment the lines below to run some tests:
# test_smiles = [
#     # True positives (examples):
#     "C[C@@H](OP(O)(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O",  # FMN-L-threonine
#     "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C",  # lumiflavin
#     "C12=NC(NC(C1=NC=3C(N2C[C@@H]([C@@H]([C@@H](COP(=O)(O)O)O)O)O)=CC(=C(C3)C)C)=O)=O",  # FMN
#     # False positive (roseoflavin-like) should not match because the benzene ring is not dimethylated as required:
#     "CN(C)c1cc2n(C[C@H](O)[C@H](O)[C@H](O)CO)c3nc(=O)[nH]c(=O)c3nc2cc1C",
# ]
#
# for smi in test_smiles:
#     result, reason = is_flavin(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")