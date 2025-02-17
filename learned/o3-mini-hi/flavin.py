"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: A derivative of the dimethylisoalloxazine 
(7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton with a substituent on the 10 position (i.e. flavin).

This version uses an updated SMARTS pattern:
  "c1c(C)cc2nc3c(nc(=O)[nH]c3=O)n(*)c2c1C"

Breakdown:
  - "c1c(C)cc2" and "c2c1C" require that the benzene ring contains two methyl groups (at positions analogous to 7 and 8)
  - "nc3c(nc(=O)[nH]c3=O)" ensures that the fused heterocycle has the two keto groups
  - "n(*)" means that the bridging nitrogen (N10) carries any substituent (not just hydrogen)
This pattern is somewhat relaxed compared to the previous version to better capture the variability in valid flavin derivatives.
"""

from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin derivative based on its SMILES string.
    A flavin is defined as a derivative of the dimethylisoalloxazine 
    (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton with a substituent at the 10 position.
    
    The classification is done in two steps:
      1. The molecule is parsed from its SMILES.
      2. An updated SMARTS pattern is used to see if the isoalloxazine core is found:
          "c1c(C)cc2nc3c(nc(=O)[nH]c3=O)n(*)c2c1C"
         This pattern requires:
           • The benzene ring shows two methyl groups.
           • The fused heterocycle has two keto groups.
           • The bridging nitrogen (position 10) carries a substituent (indicated by n(*) ).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a flavin derivative, False otherwise.
        str: Reason for the classification result.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an updated SMARTS pattern for the flavin isoalloxazine core.
    # This pattern is designed to be a bit more flexible than the previous attempt.
    flavin_smarts = "c1c(C)cc2nc3c(nc(=O)[nH]c3=O)n(*)c2c1C"
    core = Chem.MolFromSmarts(flavin_smarts)
    if core is None:
        return False, "Error constructing the flavin SMARTS pattern"
    
    # Check if the molecule contains the isoalloxazine core with the required features.
    if not mol.HasSubstructMatch(core):
        return False, "Does not contain the required dimethylisoalloxazine core with a substituent at N10"
    
    return True, "Contains isoalloxazine core (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) with a substituent at N10, classified as flavin"

# Example testing:
if __name__ == "__main__":
    test_smiles = [
        "C[C@@H](OP(O)(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O",   # FMN-L-threonine
        "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C",                                               # lumiflavin
        "Cc1cc2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3N3CCC(C)(C)c(c1C)c23",           # prenyl-FMNH2
        "O=C1N(C2=NC=3C(=CC(CO)=C(C3)C)N=C2C(N1)=O)C",                                          # 1-methyl-11-hydroxylumichrome
        "O1[C@H](N2C3=NC4=C(N=C3C(=O)NC2=O)C=C(C(=C4)C)C)[C@H](O)[C@H](O)[C@H]1CO",               # 1-(a-ribofuranosyl)-lumichrome
    ]
    for smi in test_smiles:
        result, reason = is_flavin(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")