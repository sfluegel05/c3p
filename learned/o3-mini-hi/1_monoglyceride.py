"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: CHEBI: 1-monoglyceride
Definition:
  A 1-monoglyceride is a glycerol ester where one of the two primary positions
  (i.e. one of the terminal CH2’s) on the glycerol backbone (HO–CH2–CH(OH)–CH2OH)
  is esterified via an –O–C(=O)–R bond. In this implementation two SMARTS patterns
  are used to detect the glycerol backbone with the ester function at either terminal
  position. For a valid match, the acyl substituent (mapped as 'acyl') must be aliphatic.
  
Our previous attempt failed to match many examples due to a strict SMARTS pattern.
We now relax the pattern by (a) explicitly mapping the key atoms, (b) using useChirality=False,
and (c) embedding only the necessary pattern for a glycerol backbone with one ester.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride.
    A 1-monoglyceride is defined as a glycerol ester in which the acyl chain is attached
    to one of the primary positions of the glycerol (the terminal CH2 groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a 1-monoglyceride, False otherwise.
        str: Explanation of the decision.
    """
    # Parse the input SMILES using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude molecules that contain phosphorus (uncommon in monoglycerides).
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; unlikely to be a monoglyceride"

    # Check a minimal molecular weight (heuristic).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for a typical monoglyceride"

    # Define two SMARTS patterns
    # Pattern 1: Ester group at the first terminal CH2 position.
    #   [CH2:1](O[CX3](=O)[C:acyl]) represents the first CH2 substituted with an ester oxygen,
    #   followed by a CH (with OH) and a CH2 (with free OH).
    pattern1_smarts = "[CH2:1](O[CX3](=O)[C:acyl])[CH:2](O)[CH2:3](O)"
    # Pattern 2: Ester group at the third terminal CH2 position.
    pattern2_smarts = "[CH2:1](O)[CH:2](O)[CH2:3](O[CX3](=O)[C:acyl])"
    
    pattern1 = Chem.MolFromSmarts(pattern1_smarts)
    pattern2 = Chem.MolFromSmarts(pattern2_smarts)
    if (pattern1 is None) or (pattern2 is None):
        return False, "Error in SMARTS pattern definition"

    # Helper: Build a mapping from SMARTS atom map numbers (as strings) to atom indices in the molecule.
    def get_match_mapping(query, match_tuple):
        mapping = {}
        for atom in query.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                map_num = atom.GetProp("molAtomMapNumber")
                mapping[map_num] = match_tuple[atom.GetIdx()]
        return mapping

    # Helper: Ensure that the glycerol backbone atoms (mapped as "1", "2" and "3") are acyclic.
    def backbone_is_acyclic(query, match_tuple):
        mapping = get_match_mapping(query, match_tuple)
        for label in ["1", "2", "3"]:
            if label in mapping:
                if mol.GetAtomWithIdx(mapping[label]).IsInRing():
                    return False
        return True

    # Helper: Check that the acyl substituent (mapped as "acyl") is a carbon and is not aromatic.
    def acyl_is_valid(query, match_tuple):
        mapping = get_match_mapping(query, match_tuple)
        if "acyl" not in mapping:
            return False
        acyl_atom = mol.GetAtomWithIdx(mapping["acyl"])
        if acyl_atom.GetAtomicNum() != 6:
            return False
        if acyl_atom.GetIsAromatic():
            return False
        return True

    valid_match = False
    explanation = ""

    # Use useChirality=False to relax stereochemical constraints.
    matches1 = mol.GetSubstructMatches(pattern1, useChirality=False)
    for match in matches1:
        if backbone_is_acyclic(pattern1, match) and acyl_is_valid(pattern1, match):
            valid_match = True
            explanation = "Match found: ester at terminal CH2 (position 1) of the glycerol backbone."
            break

    if not valid_match:
        matches2 = mol.GetSubstructMatches(pattern2, useChirality=False)
        for match in matches2:
            if backbone_is_acyclic(pattern2, match) and acyl_is_valid(pattern2, match):
                valid_match = True
                explanation = "Match found: ester at the other terminal CH2 (position 1) of the glycerol backbone."
                break

    if not valid_match:
        return False, "No valid glycerol backbone with ester at a primary position found or acyl substituent is aromatic."

    return True, explanation

# Optional self-test when running the script
if __name__ == "__main__":
    # A short list of examples drawn from the provided test cases.
    test_examples = [
        ("O=C(OCC(O)CO)CCCCCCCCCCCCCCC(CC)C", "AKD-2B2"),
        ("OCC(COC(CCCCCCCCC/C=C\\CCCCCCCC)=O)O", "1-(11Z-icosenoyl)glycerol"),
        ("CCCCCCCC(=O)OC[C@H](O)CO", "3-octanoyl-sn-glycerol"),
        ("CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](O)CO", "3-arachidonoyl-sn-glycerol"),
        ("C(CCCCCCCCCCCC(OCC(CO)O)=O)CCCCCCC", "1-icosanoylglycerol"),
        ("O(CC(O)CO)C(=O)C(C)=CC", "Glyceryl methylmethacrylate"),
        ("O=C(OC[C@@H](O)CO)CC1=CC=CC=C1", "1-O-(phenylacetyl)glycerol"),
    ]
    
    for smi, name in test_examples:
        result, reason = is_1_monoglyceride(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}, Reason: {reason}\n")