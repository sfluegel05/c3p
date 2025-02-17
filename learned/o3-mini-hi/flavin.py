"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: A derivative of the dimethylisoalloxazine 
(7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton 
with a substituent on the 10 position (i.e. flavin).

This version uses a refined SMARTS pattern that requires:
  - An aromatic benzene ring fused to the heterocyclic system with two explicit methyls,
  - Two keto groups on the heterocycle,
  - And the bridging (N10) nitrogen carries any non‐hydrogen substituent,
as enforced by the [!H] specification.
"""

from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin derivative based on its SMILES string.
    A flavin is defined as a derivative of the dimethylisoalloxazine
    (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton with a substituent at the 10 position.
    We enforce that:
      - The isoalloxazine core is present,
      - The benzene ring shows two methyl groups (at positions analogous to 7 and 8),
      - And the bridging nitrogen (N10) has a substituent that is not hydrogen.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a flavin derivative, False otherwise.
        str: Reason explaining the classification result.
    """
    # Parse the input SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Construct a SMARTS pattern for the isoalloxazine core:
    # Breakdown of the pattern "c1c(C)cc2nc3c(nc(=O)[nH]c3=O)n([!H])c2c1C":
    #   - "c1c(C)cc2": Starts with a benzene ring (ring opening) where one of the carbons (position 7)
    #       explicitly bears a methyl group.
    #   - "nc3c(nc(=O)[nH]c3=O)": Fused heterocyclic system containing the two keto groups.
    #   - "n([!H])": The bridging nitrogen (N10) must carry a substituent that is not hydrogen.
    #   - "c2c1C": The benzene ring is closed with another explicit methyl (position 8).
    isoalloxazine_smarts = "c1c(C)cc2nc3c(nc(=O)[nH]c3=O)n([!H])c2c1C"
    core = Chem.MolFromSmarts(isoalloxazine_smarts)
    if core is None:
        return False, "Error constructing isoalloxazine SMARTS pattern"
    
    # Check if the molecule contains the flavin isoalloxazine core.
    if not mol.HasSubstructMatch(core):
        return False, "Does not contain the required dimethylisoalloxazine core with a substituent at N10"
    
    return True, "Contains isoalloxazine core (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) with a substituent at N10, classified as flavin"

# Example testing code: (uncomment the block below to run tests)
# test_smiles = [
#     "C[C@@H](OP(O)(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O",  # FMN-L-threonine
#     "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C",  # lumiflavin
#     "CC1=C(C)C=C2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3=NC2=C1",  # FMN (illustrative)
# ]
#
# for smi in test_smiles:
#     result, reason = is_flavin(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")