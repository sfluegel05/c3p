"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: Phosphoinositide
Definition: Any phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.
That is, the molecule must contain a glycerophosphatidylinositol “core” (with two acyl chains
attached via ester bonds and a phosphate linking to an inositol ring) where at least one hydroxyl on the inositol has been phosphorylated.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    
    We use a number of heuristic checks:
      1. The molecule must be parseable.
      2. It must contain at least two phosphorus atoms (one for the glycerol linkage and at least one extra phosphoryl substituent on inositol).
      3. It must have at least two acyl ester groups (i.e. -C(=O)O- fragments).
      4. It should contain an inositol ring. Here we define an inositol ring as a six-membered aliphatic ring
         in which every ring carbon is substituted with either a hydroxyl group (OH) or a phosphate ester (OP(=O)(O)O).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check phosphorus count.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 2:
        return False, f"Only {p_count} phosphorus atom(s) found; a phosphoinositide needs at least 2 (one bridging PI and one extra phosphorylation on inositol)"
    
    # 2. Check for presence of at least two acyl ester groups.
    # This simple SMARTS pattern looks for -C(=O)O- fragments.
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Only {len(ester_matches)} acyl ester group(s) found; need at least two to represent acyl chains"
    
    # 3. Look for an inositol ring.
    # A typical myo-inositol ring is a cyclohexane with substituents. In a phosphoinositide,
    # one or more substituent positions may be phosphorylated.
    # The SMARTS below defines a six-membered ring (R6) where each carbon (CX4) is substituted with
    # either an -OH ([OX2H]) or a phosphorylated group (OP(=O)(O)O). The comma inside the square brackets
    # indicates an "OR" condition.
    inositol_smarts = ("[$(C1("
                      "[OX2H,OP(=O)(O)O])"
                      "C([OX2H,OP(=O)(O)O])"
                      "C([OX2H,OP(=O)(O)O])"
                      "C([OX2H,OP(=O)(O)O])"
                      "C([OX2H,OP(=O)(O)O])"
                      "C([OX2H,OP(=O)(O)O])1]")
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if inositol_pattern is None:
        return False, "Error generating inositol SMARTS pattern"
    
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring (with substituted O or phosphate groups) detected"
    
    # If all tests pass, then we assume the molecule is a phosphoinositide.
    return True, "Molecule has at least two phosphorus atoms, two acyl chains, and an inositol ring with phosphorylated substituents"

# Debug/test code (can be used when running as a standalone module)
if __name__ == "__main__":
    test_smiles = [
        # One of the provided examples: 1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,4-bisphosphate)
        "CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCC"
    ]
    
    for smi in test_smiles:
        result, reason = is_phosphoinositide(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")