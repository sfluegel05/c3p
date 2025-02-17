"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: isoflavones – any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton 
and its substituted derivatives.

The strategy is to use a SMARTS query that captures the critical core:
    c1ccc2c(c1)oc(=O)c(c2)c3ccccc3

This pattern describes a benzopyran-4-one (chromen-4-one) where the pyran (ring C) bears an aryl substituent at its 3‐position.
Substructure matching via RDKit then allows extra substituents and decorations (e.g. methoxy groups, glycosidic linkages),
so that many substituted isoflavones are still detected.
"""

from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone (3-aryl-1-benzopyran-4-one) based on its SMILES string.
    
    The approach uses a SMARTS pattern designed to match the isoflavone core.
    The pattern used is:
         c1ccc2c(c1)oc(=O)c(c2)c3ccccc3
    which searches for:
      • A benzopyran-4-one (chromen-4-one) framework (the fused benzene + pyran-one)
      • With a phenyl (aryl) substituent attached to the pyran ring (at its 3‑position)
    
    Args:
         smiles (str): SMILES string of the molecule.
         
    Returns:
         bool: True if the molecule contains an isoflavone skeleton; False otherwise.
         str: A reason explaining the classification decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the isoflavone core.
    # This pattern looks for a fused benzene-pyranone ring system with an attached aryl group.
    isoflavone_smarts = "c1ccc2c(c1)oc(=O)c(c2)c3ccccc3"
    core = Chem.MolFromSmarts(isoflavone_smarts)
    if core is None:
        return False, "Failed to compile the SMARTS pattern"
    
    # Check if the molecule contains the isoflavone core.
    if mol.HasSubstructMatch(core):
        return True, "Molecule contains a 3-aryl-1-benzopyran-4-one skeleton and its substituted derivatives."
    else:
        return False, "Molecule does not contain a suitable 3-aryl-1-benzopyran-4-one skeleton."

# Example test cases (a few examples from the list provided)
if __name__ == "__main__":
    test_cases = {
        "formononetin": "COc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",
        "Lupiwighteone": "O1C2=C(CC=C(C)C)C(O)=CC(O)=C2C(=O)C(C3=CC=C(O)C=C3)=C1",
        "2',4',5,7-Tetrahydroxy-8-prenylisoflavone": "O1C2=C(CC=C(C)C)C(O)=CC(O)=C2C(=O)C(C3=C(O)C=C(O)C=C3)=C1",
        "5,6,7-trimethoxy-3-(3,4,5-trimethoxyphenyl)-1-benzopyran-4-one": "COC1=CC(=CC(=C1OC)OC)C2=COC3=CC(=C(C(=C3C2=O)OC)OC)OC",
        "2'-hydroxydaidzein": "Oc1ccc(c(O)c1)-c1coc2cc(O)ccc2c1=O",
        "glycitein": "COc1cc2c(cc1O)occ(-c1ccc(O)cc1)c2=O",
        "millewanin H": "CC(C)=CCc1c(O)c(CC(O)C(C)=C)c2occ(-c3ccc(O)c(O)c3)c(=O)c2c1O",
        "7,4'-dihydroxy-3'-methoxyisoflavone": "COc1cc(ccc1O)-c1coc2cc(O)ccc2c1=O",
        "bowdichione": "COC1=CC(=O)C(=CC1=O)c1coc2cc(O)ccc2c1=O",
        # A false positive example (should return False): 
        "quercetin 3,3'-bissulfate": "Oc1ccc(cc1OS([O-])(=O)=O)-c1oc2cc([O-])cc(O)c2c(=O)c1OS([O-])(=O)=O",
    }
    for name, smi in test_cases.items():
        result, reason = is_isoflavones(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}, Reason: {reason}\n")