"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: CHEBI: 1-monoglyceride 
Definition: A monoglyceride in which the acyl substituent is located at a primary position (i.e. at position 1 or 3) on glycerol.
Expected features:
    - Exactly one ester group (O–C(=O)–R)
    - A glycerol backbone where one of the terminal CH2 groups is esterified while the other two hydroxyl groups remain free.
Note: This approach uses substructure SMARTS searches and some heuristic checks.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if the given SMILES string represents a 1-monoglyceride,
    i.e. a monoglyceride with the acyl substituent attached at a primary position (sn-1 or sn-3)
    of the glycerol backbone.
    
    We use the following criteria:
    1. The molecule must be parsed successfully.
    2. The molecule must contain exactly one ester group (pattern: [CX3](=O)O).
    3. The molecule must contain a glycerol backbone pattern in which one terminal CH2 is esterified 
       and the other two hydroxyl groups are available.
    
    We search for two possible patterns (the acyl group may be on either end):
      Pattern 1:  [CH2](OC(=O)*)[CH](O)[CH2](O)
      Pattern 2:  [CH2](O)[CH](O)[CH2](OC(=O)*)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 1-monoglyceride, False otherwise.
        str: A message providing the reason for the classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an ester SMARTS pattern (acyl carbon and the ester oxygen)
    ester_smarts = "[CX3](=O)O"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected exactly one ester group, found {len(ester_matches)}"
    
    # Define two possible glycerol backbone patterns where one terminal CH2 is esterified.
    # Pattern 1: acyl group at the first position
    pattern1_smarts = "[CH2](OC(=O)*)[CH](O)[CH2](O)"
    pattern1 = Chem.MolFromSmarts(pattern1_smarts)
    # Pattern 2: acyl group at the third (other terminal) position
    pattern2_smarts = "[CH2](O)[CH](O)[CH2](OC(=O)*)"
    pattern2 = Chem.MolFromSmarts(pattern2_smarts)
    
    match1 = mol.HasSubstructMatch(pattern1)
    match2 = mol.HasSubstructMatch(pattern2)
    
    if not (match1 or match2):
        return False, "Glycerol backbone with ester at a primary position not found"
    
    # Optionally, we can perform additional sanity checks.
    # For example, we expect a certain size in molecular weight in typical monoglycerides.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for a typical monoglyceride"
    
    return True, "Found a glycerol backbone with one ester group at a primary position"

# For testing purposes, you can run a few examples:
if __name__ == "__main__":
    test_smiles = [
        # 1-(11Z-icosenoyl)glycerol
        "OCC(COC(CCCCCCCCC/C=C\\CCCCCCCC)=O)O",
        # 3-octanoyl-sn-glycerol
        "CCCCCCCC(=O)OC[C@H](O)CO",
        # rac-1-monolauroylglycerol
        "CCCCCCCCCCCC(=O)OCC(O)CO",
        # 1-stearoyl-sn-glycerol
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)CO",
        # Example of non-monoglyceride
        "CC(=O)OC1=CC=CC=C1"
    ]
    for s in test_smiles:
        valid, reason = is_1_monoglyceride(s)
        print(f"SMILES: {s}\n Classified: {valid}, Reason: {reason}\n")