"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: CHEBI: 1-monoglyceride
Definition: A monoglyceride in which the acyl substituent is located at a primary (non-cyclic) position on a glycerol backbone.
This implementation uses two SMARTS patterns to detect a glycerol backbone (HO–CH2–CH(OH)–CH2OH)
in which one of the terminal CH2 groups is esterified (i.e. linked via –O–C(=O)–R).
The acyl substituent (labeled "AC") is then checked to be an aliphatic (non-aromatic) carbon.
A few additional checks (e.g. minimal molecular weight, absence of phosphorus) are used to exclude unlikely cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride.
    A 1-monoglyceride is defined as a glycerol ester in which the acyl chain is attached
    to one of the primary positions of the glycerol (i.e. at CH2)
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a 1-monoglyceride, False otherwise.
      str: Explanation of the decision.
    """
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules with phosphorus (unlikely in monoglycerides)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely not a monoglyceride"
        
    # Minimal molecular weight check (heuristic)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for a typical monoglyceride"
    
    # Define two SMARTS for glycerol backbone with an ester linkage at a terminal CH2.
    # In the SMARTS:
    #   - [CH2:1] and [CH2:3] are the two terminal carbons.
    #   - [CH:2] is the central carbon.
    #   - In pattern1, the ester is on the first CH2:  [CH2:1](O[CX3](=O)[C:AC])
    #   - In pattern2, the ester is on the third CH2:  [CH2:3](O[CX3](=O)[C:AC])
    #
    # Pattern 1: Ester at the first CH2 group:
    pattern1_smarts = "[CH2:1](O[CX3](=O)[C:AC])[CH:2](O)[CH2:3](O)"
    # Pattern 2: Ester at the third CH2 group:
    pattern2_smarts = "[CH2:1](O)[CH:2](O)[CH2:3](O[CX3](=O)[C:AC])"
    
    pattern1 = Chem.MolFromSmarts(pattern1_smarts)
    pattern2 = Chem.MolFromSmarts(pattern2_smarts)
    
    if (pattern1 is None) or (pattern2 is None):
        return False, "Error in SMARTS pattern definition"
    
    # Helper: get dictionary mapping SMARTS map labels to the corresponding mol atom indices
    def get_match_mapping(query, match_tuple):
        mapping = {}
        for atom in query.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                map_num = atom.GetProp("molAtomMapNumber")
                mapping[map_num] = match_tuple[atom.GetIdx()]
        return mapping

    # Helper: make sure the glycerol backbone atoms (labeled "1", "2", "3") are not in a ring
    def backbone_is_acyclic(query, match_tuple):
        mapping = get_match_mapping(query, match_tuple)
        for label in ["1", "2", "3"]:
            if label in mapping:
                if mol.GetAtomWithIdx(mapping[label]).IsInRing():
                    return False
        return True

    # Helper: check that the acyl substituent (mapped as "AC") is an aliphatic carbon (i.e. not aromatic)
    def acyl_is_valid(query, match_tuple):
        mapping = get_match_mapping(query, match_tuple)
        if "AC" not in mapping:
            return False
        acyl_atom = mol.GetAtomWithIdx(mapping["AC"])
        if acyl_atom.GetAtomicNum() != 6:
            return False
        if acyl_atom.GetIsAromatic():
            return False
        return True

    valid_match = False
    explanation = ""
    
    # Try matching using pattern 1 (ester at first CH2)
    matches1 = mol.GetSubstructMatches(pattern1)
    for match in matches1:
        if backbone_is_acyclic(pattern1, match) and acyl_is_valid(pattern1, match):
            valid_match = True
            explanation = "Match found with ester at the first CH2 position (primary position)."
            break

    # If pattern 1 did not produce a valid match, try pattern 2 (ester at third CH2)
    if not valid_match:
        matches2 = mol.GetSubstructMatches(pattern2)
        for match in matches2:
            if backbone_is_acyclic(pattern2, match) and acyl_is_valid(pattern2, match):
                valid_match = True
                explanation = "Match found with ester at the third CH2 position (primary position)."
                break

    if not valid_match:
        return False, "No valid glycerol backbone with ester at a primary position found or acyl substituent is aromatic."
    
    return True, explanation


# Optional: simple tests when running the script
if __name__ == "__main__":
    test_examples = [
        ("O=C(OCC(O)CO)CCCCCCCCCCCCCCC(CC)C", "AKD-2B2"),
        ("OCC(COC(CCCCCCCCC/C=C\\CCCCCCCC)=O)O", "1-(11Z-icosenoyl)glycerol"),
        ("CCCCCCCC(=O)OC[C@H](O)CO", "3-octanoyl-sn-glycerol"),
        ("CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](O)CO", "3-arachidonoyl-sn-glycerol"),
        ("C(CCCCCCCCCCCC(OCC(CO)O)=O)CCCCCCC", "1-icosanoylglycerol"),
        ("O(CC(O)CO)C(=O)C(C)=CC", "Glyceryl methylmethacrylate"),
        ("CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](O)CO", "3-oleoyl-sn-glycerol"),
        ("O(C(=O)CCCCCCCCCCCCC/C=C\\CCCCCCCC)C[C@@H](O)CO", "MG(24:1(15Z)/0:0/0:0)"),
        ("[H][C@@]12CC[C@]3(CC1=C)C=CC(=O)[C@@](C)(CCC(=O)OCC(O)CO)[C@]3([H])C2", "Platencinic acid glycerol ester"),
        ("CCCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@@H](O)CO", "1-linoleoyl-sn-glycerol"),
        ("O=C(OC[C@@H](O)CO)CC1=CC=CC=C1", "1-O-(phenylacetyl)glycerol"),
        ("CCCCCCCCCCC(=O)OC[C@@H](O)CO", "rac-1-monolauroylglycerol"),
        ("CCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)CO", "1-stearoyl-sn-glycerol"),
    ]
    
    for smi, name in test_examples:
        result, reason = is_1_monoglyceride(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}, Reason: {reason}\n")