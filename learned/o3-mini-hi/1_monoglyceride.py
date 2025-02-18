"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: CHEBI: 1-monoglyceride
Definition: A monoglyceride in which the acyl substituent is located at a primary position (i.e. at position 1 or equivalently 3)
This improved version:
  1. Parses the SMILES string.
  2. Excludes molecules that contain phosphorus.
  3. Requires exactly one ester group (using the SMARTS "[CX3](=O)O").
  4. Searches for a glycerol backbone with the ester at a primary (terminal) CH2;
     we define two patterns – one with the acyl substitution at the first carbon and one with it at the third.
     In these SMARTS the backbone carbons are explicitly labeled (using :1, :2, :3).
  5. Checks that the three backbone atoms are not part of any ring (to reduce false positives).
  6. Requires a minimal molecular weight (100 Da).
This heuristic approach improves specificity by verifying that the glycerol backbone appears as a freely hanging (acyclic) unit.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if the given SMILES string represents a 1-monoglyceride
    in which the acyl substituent is located at a primary position of the glycerol backbone.
    
    Args:
      smiles (str): SMILES representation of the molecule.
      
    Returns:
      bool: True if the molecule is classified as a 1-monoglyceride
      str: An explanation for the decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules containing phosphorus (atomic number 15)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely not a monoglyceride"
        
    # Count ester groups using the SMARTS pattern "[CX3](=O)O"
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected exactly one ester group, found {len(ester_matches)}"
    
    # Define two glycerol backbone patterns with explicit atom mapping.
    # In pattern1 the acyl group (i.e. the ester) is on the first CH2 (position 1).
    # In pattern2 the ester is on the last CH2 (position 3), which is equivalent.
    #
    # The mapping labels (:1, :2, :3) mark the three backbone carbons.
    pattern1_smarts = "[CH2:1](OC(=O)[*])[CH:2](O)[CH2:3](O)"
    pattern2_smarts = "[CH2:1](O)[CH:2](O)[CH2:3](OC(=O)[*])"
    pattern1 = Chem.MolFromSmarts(pattern1_smarts)
    pattern2 = Chem.MolFromSmarts(pattern2_smarts)
    
    # Get matches for both glycerol backbone patterns.
    matches1 = mol.GetSubstructMatches(pattern1)
    matches2 = mol.GetSubstructMatches(pattern2)
    
    # A helper function to check that the backbone carbons (mapped as "1", "2" and "3")
    # are not in a ring. This ensures that the glycerol backbone is not incorporated in a cyclic system,
    # which was a common cause of false positives.
    def backbone_is_acyclic(query, match_tuple):
        # Build a mapping from the query atom mapping number to the corresponding atom in mol.
        mapping = {}
        for atom in query.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                map_num = atom.GetProp("molAtomMapNumber")
                mapping[map_num] = match_tuple[atom.GetIdx()]
        # Check that each of backbone carbons labeled "1", "2" and "3" is not in a ring.
        for label in ["1", "2", "3"]:
            if label in mapping:
                if mol.GetAtomWithIdx(mapping[label]).IsInRing():
                    return False
        return True
    
    # Check if any match from pattern1 or pattern2 has an acyclic backbone.
    valid_backbone = False
    for match in matches1:
        if backbone_is_acyclic(pattern1, match):
            valid_backbone = True
            break
    if not valid_backbone:
        for match in matches2:
            if backbone_is_acyclic(pattern2, match):
                valid_backbone = True
                break

    if not valid_backbone:
        return False, "Glycerol backbone with ester at a primary position not found or is part of a cyclic structure"
    
    # Check minimal molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for a typical monoglyceride"
    
    return True, "Found glycerol backbone with one ester group at a primary (non-cyclic) position"

# Optional: Simple testing when run as script
if __name__ == "__main__":
    test_examples = [
        # True positives:
        "O=C(OCC(O)CO)CCCCCCCCCCCCCCC(CC)C",                # AKD-2B2
        "OCC(COC(CCCCCCCCC/C=C\\CCCCCCCC)=O)O",               # 1-(11Z-icosenoyl)glycerol
        "CCCCCCCC(=O)OC[C@H](O)CO",                          # 3-octanoyl-sn-glycerol
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](O)CO",            # 3-arachidonoyl-sn-glycerol
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)CO",               # 1-stearoyl-sn-glycerol
        "O(C[C@@H](O)CO)C(=O)C",                             # (R)-glycerol 1-acetate
        # False positive examples (should be rejected):
        "O=C(OC[C@@H](O)CO)C1=C(O)C(=CC=C1C)",               # Glycerol 1-hydroxy-2,5-dimethyl benzoate
    ]
    for s in test_examples:
        result, explanation = is_1_monoglyceride(s)
        print(f"SMILES: {s}\n Classified: {result}, Reason: {explanation}\n")