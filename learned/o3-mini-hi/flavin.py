"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: A derivative of the dimethylisoalloxazine 
(7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton with a substituent at the 10 position (i.e. flavin).

This version uses an updated SMARTS pattern that maps the bridging N (position 10) and requires it to have a non‐hydrogen substituent.
The SMARTS pattern used is:
  "c1c(C)cc2nc3c(nc(=O)[nH]c3=O)[n:1]([!H])c2c1C"
Breakdown:
  - "c1c(C)cc" and then later "c2c1C" enforce that the benzene ring has two methyl groups (positions 7 and 8).
  - "nc3c(nc(=O)[nH]c3=O)" captures the fused heterocycle with the two keto functions.
  - "[n:1]([!H])" explicitly maps the central (N10) nitrogen and requires that it carries at least one substituent that is not hydrogen.
"""

from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is classified as a flavin derivative based on its SMILES string.
    A flavin derivative is defined as a derivative of the dimethylisoalloxazine 
    (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton with a substituent at the 10 position.
    
    The procedure is:
      1. Parse the molecule from the SMILES.
      2. Apply a SMARTS substructure search with an updated pattern designed to capture the isoalloxazine core
         and specifically require that the bridging nitrogen (N10) carries a substituent (i.e. not just hydrogen).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as flavin.
        str: Explanation for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an updated SMARTS pattern for the isoalloxazine core.
    # Here we force:
    #  • a benzene ring with two methyl groups (positions analogous to 7 and 8),
    #  • a fused heterocycle having two keto groups,
    #  • and the bridging nitrogen (N10) [mapped as [n:1]] which carries a non-hydrogen substituent ([!H]).
    flavin_smarts = "c1c(C)cc2nc3c(nc(=O)[nH]c3=O)[n:1]([!H])c2c1C"
    core = Chem.MolFromSmarts(flavin_smarts)
    if core is None:
        return False, "Error constructing the flavin SMARTS pattern"
    
    # Check if the molecule contains the flavin isoalloxazine core with the required N10 substitution.
    if not mol.HasSubstructMatch(core):
        return False, "Does not contain the required dimethylisoalloxazine core with a substituent at N10"
    
    return True, "Contains isoalloxazine core (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) with substituent at N10, classified as flavin"

# Example test cases (can be removed or commented out when using this code as a module):
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