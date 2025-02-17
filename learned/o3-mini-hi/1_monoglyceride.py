"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: CHEBI: 1-monoglyceride
Definition: A monoglyceride in which the acyl substituent is located at a primary position (i.e. at position 1 or, equivalently, position 3) on glycerol.
Criteria used:
  1. The SMILES must be parsed successfully.
  2. The molecule must NOT contain phosphorus (to avoid phospholipids).
  3. The molecule must have exactly one ester group, as identified by the SMARTS pattern “[CX3](=O)O”.
  4. A glycerol backbone must be detected with one of two patterns:
         Pattern 1: acyl group at the first (sn-1) position:
             [CH2](OC(=O)[*])[CH](O)[CH2](O)
         Pattern 2: acyl group at the third (sn-3) position:
             [CH2](O)[CH](O)[CH2](OC(=O)[*])
  5. A minimal molecular weight is required (set to 100 Da) so that valid, short-chain derivatives are not rejected.
  
Note: This approach is heuristic and relies on substructure SMARTS searches.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if the given SMILES string represents a 1-monoglyceride,
    i.e. a monoglyceride where the acyl substituent is at a primary (sn-1 or sn-3) position.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if it is likely a 1-monoglyceride, False otherwise.
        str: Explanation of the decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules containing phosphorus (P) to reduce false positives (phospholipids, etc.)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely not a monoglyceride"

    # Count ester groups using the SMARTS "[CX3](=O)O"
    ester_smarts = "[CX3](=O)O"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected exactly one ester group, found {len(ester_matches)}"
    
    # Define two glycerol backbone patterns where one terminal CH2 is esterified.
    # The '*' in the pattern means any R group attached to the carbonyl oxygen.
    pattern1_smarts = "[CH2](OC(=O)[*])[CH](O)[CH2](O)"
    pattern2_smarts = "[CH2](O)[CH](O)[CH2](OC(=O)[*])"
    pattern1 = Chem.MolFromSmarts(pattern1_smarts)
    pattern2 = Chem.MolFromSmarts(pattern2_smarts)
    
    match1 = mol.HasSubstructMatch(pattern1)
    match2 = mol.HasSubstructMatch(pattern2)
    if not (match1 or match2):
        return False, "Glycerol backbone with ester at a primary position not found"
    
    # Check minimal molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for a typical monoglyceride"
    
    return True, "Found glycerol backbone with one ester group at a primary position"

# Testing examples (optional)
if __name__ == "__main__":
    # A list of example SMILES strings for molecules that should be classified as 1-monoglycerides.
    test_smiles = [
        "O=C(OCC(O)CO)CCCCCCCCCCCCCCC(CC)C",  # AKD-2B2
        "OCC(COC(CCCCCCCCC/C=C\\CCCCCCCC)=O)O",  # 1-(11Z-icosenoyl)glycerol
        "CCCCCCCC(=O)OC[C@H](O)CO",             # 3-octanoyl-sn-glycerol
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](O)CO",# 3-arachidonoyl-sn-glycerol
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)CO",   # 1-stearoyl-sn-glycerol
        "O(C[C@@H](O)CO)C(=O)C",                 # (R)-glycerol 1-acetate (should pass despite low weight)
        "O(CC(O)CO)C(=O)CC"                     # Glycerol 1-propanoate (will fail if molecular weight is below threshold)
    ]
    for s in test_smiles:
        valid, reason = is_1_monoglyceride(s)
        print(f"SMILES: {s}\n Classified: {valid}, Reason: {reason}\n")