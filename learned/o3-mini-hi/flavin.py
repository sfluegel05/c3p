"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: A derivative of the dimethylisoalloxazine 
(7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton with a substituent 
at the 10 position (i.e. flavin).

This version uses a revised SMARTS pattern. In our updated pattern we require:
  • a fused bicyclic/heterocyclic system with a benzene ring substituted with two methyl groups (7,8-dimethyl),
  • the presence of two keto groups on the fused heterocycle,
  • and that the bridging nitrogen (at the 10 position) carries a substituent (i.e. is not just –H).

The final SMARTS pattern is:
  "c1c(C)cc2nc3c(nc(=O)[nH]c3=O)[n]([!H])c2c1C"
This pattern is somewhat relaxed with respect to atom mapping to allow some variation in how the core is drawn.
"""
from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is classified as a flavin derivative based on its SMILES string.
    A flavin derivative is defined as a derivative of the dimethylisoalloxazine 
    (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton with a substituent at the 10 position.
    
    The procedure is:
      1. Parse the molecule from the SMILES.
      2. Use a SMARTS substructure search with our updated pattern designed to capture the isoalloxazine core.
         This pattern enforces:
           - A benzene ring with two methyl groups (mimicking substitutions at positions 7 and 8),
           - a fused heterocycle containing two keto groups,
           - and a bridging nitrogen that carries at least one non‐hydrogen substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a flavin derivative, False otherwise.
        str: Explanation for the classification.
    """
    
    # Parse the SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Updated SMARTS for the dimethylisoalloxazine core with substituent at N10.
    # Explanation of the pattern:
    #  • "c1c(C)cc" and "c2c1C" enforces a benzene-like ring bearing two methyl groups.
    #  • "nc3c(nc(=O)[nH]c3=O)" captures the fused heterocycle with two keto groups.
    #  • "[n]([!H])" insists that the bridging nitrogen has a substituent that is not hydrogen.
    flavin_smarts = "c1c(C)cc2nc3c(nc(=O)[nH]c3=O)[n]([!H])c2c1C"
    
    core = Chem.MolFromSmarts(flavin_smarts)
    if core is None:
        return False, "Error constructing the flavin SMARTS pattern"
    
    # Check if the molecule contains the isoalloxazine (flavin) core.
    if not mol.HasSubstructMatch(core):
        return False, "Does not contain the required dimethylisoalloxazine core with a substituent at N10"
    
    return True, "Contains isoalloxazine core with two methyl groups and N10 substitution, classified as flavin"


# Example test cases to illustrate usage (can be removed/commented out if using as a module):
if __name__ == "__main__":
    test_smiles = [
        "C[C@@H](OP(O)(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O",  # FMN-L-threonine
        "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C",                                               # lumiflavin
        "Cc1cc2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3N3CCC(C)(C)c(c1C)c23",           # prenyl-FMNH2
        "O=C1N(C2=NC=3C(=CC(CO)=C(C3)C)N=C2C(N1)=O)C",                                          # 1-methyl-11-hydroxylumichrome
        "O1[C@H](N2C3=NC4=C(N=C3C(=O)NC2=O)C=C(C(=C4)C)C)[C@H](O)[C@H](O)[C@H]1CO",               # 1-(a-ribofuranosyl)-lumichrome
        "CC1=C(C)C=C2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3=NC2=C1",                  # riboflavin
    ]
    
    for smi in test_smiles:
        result, reason = is_flavin(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")