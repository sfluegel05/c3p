"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
"""
Classifies: Phosphatidylinositol (glycerophosphoinositol)
Definition: Any glycerophosphoinositol having one phosphatidyl group esterified
to one of the hydroxy groups of inositol.
The logic here is to look for:
  1) an inositol-phosphate “headgroup” (a six‐membered ring with multiple OH groups and 
     a P–O linkage) and 
  2) evidence of a diacylglycerol segment (by having at least two ester bonds, as found in
     fatty acyl chains attached to glycerol).
Note:
Due to the variability of phosphatidylinositol structures, the SMARTS patterns used are
heuristic and may not catch every valid example.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is of the class phosphatidylinositol (glycerophosphoinositol)
    by checking if it contains an inositol-phosphate headgroup and a diacylglycerol
    (phosphatidyl) moiety as evidenced by at least two ester groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as phosphatidylinositol, False otherwise
        str: Explanation for the classification decision
    """
    # Try to parse the provided SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ------------------------------------------------------------------------------------------------
    # 1. Check for inositol-phosphate headgroup.
    # 
    # The inositol headgroup is a cyclohexane with -OH substituents on almost every carbon.
    # In most phosphatidylinositols one of the -OH groups is replaced by a phosphate ester;
    # here we look for a substructure similar to:
    #      C1(C(C(C(C(C1O)O)O)O)O)OP
    # This pattern (ignoring possible stereochemistry) is meant to capture a typical inositol
    # ring bearing a phosphate group.
    # ------------------------------------------------------------------------------------------------
    inositol_phosphate_smarts = "C1(C(C(C(C(C1O)O)O)O)O)OP"
    inositol_pattern = Chem.MolFromSmarts(inositol_phosphate_smarts)
    if inositol_pattern is None:
        return False, "Error creating inositol SMARTS pattern"
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol-phosphate headgroup found"

    # ------------------------------------------------------------------------------------------------
    # 2. Check for the diacylglycerol moiety (phosphatidyl group).
    #
    # In a phosphatidylinositol the phosphate is esterified to a glycerol segment that carries
    # two fatty acyl chains. We use a heuristic: the presence of at least two ester bonds ([O][C](=O))
    # is expected due to the two fatty acid moieties.
    # ------------------------------------------------------------------------------------------------
    ester_smarts = "[OX2][CX3](=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return False, "Error creating ester SMARTS pattern"
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Only {len(ester_matches)} ester group(s) found; need at least 2 evidencing diacyl chains"

    # Optionally, one might add further checks here (e.g. verifying a glycerol backbone
    # or the linkage from phosphate to both an inositol group and a diacylglycerol segment). 
    # For the present heuristic two requirements are used.

    return True, "Contains an inositol-phosphate headgroup with diacylglycerol-esters consistent with phosphatidylinositol"

# Example usage:
if __name__ == "__main__":
    # Test one of the provided examples:
    test_smiles = "P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCC(O)/C=C/C=O)COC(=O)CCCCCCCCCCCCCCC)(O)=O"
    valid, reason = is_phosphatidylinositol(test_smiles)
    print(valid, reason)