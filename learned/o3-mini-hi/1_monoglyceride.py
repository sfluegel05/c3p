"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: CHEBI: 1-monoglyceride
Definition: A monoglyceride in which the acyl substituent is located at a primary (non-cyclic) position on a glycerol backbone.
Improvements in this version:
  • Do not require a global count of ester groups.
  • Look specifically for a glycerol backbone (HO-CH2-CH(OH)-CH2-OH) in which one of the terminal CH2 groups has an ester linkage.
  • Use SMARTS that explicitly maps the ester linkage and the acyl substituent.
  • Check that the glycerol backbone atoms are acyclic.
  • Check that the acyl substituent is aliphatic (i.e. its first carbon, right after the carbonyl, is not aromatic).
  • Require a minimal molecular weight.
  
Note: This is a heuristic approach and may not perfectly separate all cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if the given SMILES string represents a 1-monoglyceride --
    i.e. a glycerol ester in which the acyl chain is attached to a primary position (position 1 or equivalently position 3).
    
    The approach:
      1. Parse the SMILES and reject if invalid or if phosphorus is present.
      2. Use two backbone SMARTS patterns with explicit atom mappings.
         The patterns include the glycerol backbone and the ester group. In the SMARTS the acyl substituent (after the carbonyl)
         is mapped (as "AC") so that we can check that that carbon is aliphatic.
         Pattern1 has the ester at the first CH2; pattern2 has the ester at the last CH2.
      3. For any obtained match, ensure that the three backbone atoms (labeled "1", "2", and "3") are not in a ring.
      4. For the match, retrieve the acyl substituent atom (mapping "AC") and reject if it is aromatic.
      5. Check that the overall molecular weight is above a minimal threshold.
      
    Args:
      smiles (str): SMILES representation of the molecule.
      
    Returns:
      bool: True if the molecule is judged to be a 1-monoglyceride.
      str: Explanation of the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules with phosphorus (atomic number 15)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely not a monoglyceride"
        
    # Minimal molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for a typical monoglyceride"
        
    # Define two glycerol backbone SMARTS with explicit mapping.
    # In these SMARTS:
    #   - The glycerol backbone carbons are labeled ":1", ":2", and ":3".
    #   - The ester linkage is explicitly defined.
    #   - The acyl substituent (just after the carbonyl) is labeled as "AC".
    #
    # Pattern 1: Ester at the first CH2 group.
    #   [CH2:1](O[C:OS](=[O:OX])[C:AC])[CH:2](O)[CH2:3](O)
    #
    # Pattern 2: Ester at the third CH2 group.
    #   [CH2:1](O)[CH:2](O)[CH2:3](O[C:OS](=[O:OX])[C:AC])
    #
    pattern1_smarts = "[CH2:1](O[C:OS](=[O:OX])[C:AC])[CH:2](O)[CH2:3](O)"
    pattern2_smarts = "[CH2:1](O)[CH:2](O)[CH2:3](O[C:OS](=[O:OX])[C:AC])"
    
    pattern1 = Chem.MolFromSmarts(pattern1_smarts)
    pattern2 = Chem.MolFromSmarts(pattern2_smarts)
    
    if (pattern1 is None) or (pattern2 is None):
        return False, "Error in SMARTS definition"
    
    # Helper function: Given a query and a match tuple,
    # verify that the backbone atoms labeled "1", "2" and "3" are not in a ring.
    # Also, retrieve a dictionary mapping the label (from the SMARTS 'molAtomMapNumber' property)
    # to the corresponding atom index in the molecule.
    def get_match_mapping(query, match_tuple):
        mapping = {}
        for atom in query.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                map_num = atom.GetProp("molAtomMapNumber")
                # Note: atom.GetIdx() is the index in the query; we use that to look up the corresponding mol atom index.
                mapping[map_num] = match_tuple[atom.GetIdx()]
        return mapping

    def backbone_is_acyclic(query, match_tuple):
        mapping = get_match_mapping(query, match_tuple)
        for label in ["1", "2", "3"]:
            if label in mapping:
                if mol.GetAtomWithIdx(mapping[label]).IsInRing():
                    return False
        return True
    
    # Helper function: Check that the "acyl" atom (mapped as "AC") is aliphatic.
    def acyl_is_valid(query, match_tuple):
        mapping = get_match_mapping(query, match_tuple)
        # We expect the acyl atom to have the label "AC"
        if "AC" not in mapping:
            return False
        acyl_atom = mol.GetAtomWithIdx(mapping["AC"])
        # Require that the first acyl atom is a carbon and is not aromatic.
        if acyl_atom.GetAtomicNum() != 6:
            return False
        if acyl_atom.GetIsAromatic():
            return False
        return True

    valid_backbone = False
    # Try pattern1 matches
    matches1 = mol.GetSubstructMatches(pattern1)
    for match in matches1:
        if backbone_is_acyclic(pattern1, match) and acyl_is_valid(pattern1, match):
            valid_backbone = True
            break
    # If not found, try pattern2 matches.
    if not valid_backbone:
        matches2 = mol.GetSubstructMatches(pattern2)
        for match in matches2:
            if backbone_is_acyclic(pattern2, match) and acyl_is_valid(pattern2, match):
                valid_backbone = True
                break

    if not valid_backbone:
        return False, "No valid glycerol backbone with ester at a primary (non-cyclic) position found or acyl substituent is aromatic"
    
    return True, "Found glycerol backbone with one ester group at a primary (non-cyclic) position and an acceptable acyl substituent"

# Optional: run some simple tests if executed as a script
if __name__ == "__main__":
    test_examples = [
        # True positives:
        ("O=C(OCC(O)CO)CCCCCCCCCCCCCCC(CC)C", "AKD-2B2"),
        ("OCC(COC(CCCCCCCCC/C=C\\CCCCCCCC)=O)O", "1-(11Z-icosenoyl)glycerol"),
        ("CCCCCCCC(=O)OC[C@H](O)CO", "3-octanoyl-sn-glycerol"),
        ("CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](O)CO", "3-arachidonoyl-sn-glycerol"),
        ("CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)CO", "1-stearoyl-sn-glycerol"),
        ("O(C[C@@H](O)CO)C(=O)C", "(R)-glycerol 1-acetate"),
        ("O=C(OC[C@@H](O)CO)/C=C/CCCCCCC(=O)OC", "Cladosporester C"),  # previously false negative
        ("O=C(OC[C@@H](O)CO)CC1=CC=CC=C1", "1-O-(phenylacetyl)glycerol"),
        # False positives (should be rejected):
        ("O=C(OC[C@@H](O)CO)C1=C(O)C(=CC=C1C)", "Glycerol 1-hydroxy-2,5-dimethyl benzoate"),
    ]
    
    for smi, name in test_examples:
        result, explanation = is_1_monoglyceride(smi)
        print(f"Name: {name}\nSMILES: {smi}\nClassified: {result}, Reason: {explanation}\n")