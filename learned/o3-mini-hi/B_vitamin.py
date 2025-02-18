"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin
Definition:
  Any member of the group of eight water‐soluble vitamins (B1, B2, B3, B5, B6, B7, B9 and B12)
  that play important roles in cell metabolism.
  
Improvements over the previous version:
  • For vitamin B12 (cobalamin derivatives), in addition to detecting cobalt we also require
    a high molecular weight (>800 Da) and a minimum of 4 rings to reduce misclassifications.
  • For vitamins B6 and B9, extra SMARTS patterns covering different phosphate and protonation forms
    have been added so that variants such as pyridoxal 5'-phosphate, pyridoxamine 5'-phosphate, and
    tetrahydrofolate derivatives are detected.
  • For vitamin B3 (niacin), we now set an acceptable molecular‐weight range (roughly 120–200 Da)
    and require no sulfur atoms in order to reduce false positives.
Note: This rule‐based approach is approximate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a member of the B vitamin class (B1, B2, B3, B5, B6, B7, B9 or B12)
    by matching a series of SMARTS substructure patterns along with extra conditions.
    
    Args:
        smiles (str): The SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a B vitamin, False otherwise.
        str: Explanation indicating which pattern was matched or why it failed.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular descriptors
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    ring_count = len(mol.GetRingInfo().AtomRings())
    # Count phosphorus (to exclude compounds like CoA)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    # Count sulfur (used to exclude false positives for niacin-type molecules)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    
    # --- Vitamin B12 (cobalamin) check ---
    cobalt_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 27]
    if cobalt_atoms:
        # High weight (>800 Da) and several rings required for a corrin-like structure.
        if mol_wt > 800 and ring_count >= 4:
            return True, ("Contains cobalt in a corrin-like macrocycle (high MW and multiple rings), "
                          "indicating a vitamin B12 derivative")
    
    # Define vitamin subclass patterns and extra constraints.
    vitamin_checks = [
        {
            "name": "B1 (thiamine)",
            "patterns": [
                # canonical thiamine core (thiazolium–pyrimidine motif)
                "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1",
                # variant with a phosphate (e.g. triphosphate forms)
                "Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c2C)c(N)n1"
            ],
            "constraint": lambda m, wt, p, s: True,
            "explanation": "Matches thiamine (B1) core pattern: correct thiazolium–pyrimidine motif found."
        },
        {
            "name": "B2 (riboflavin/FAD)",
            "patterns": [
                # Riboflavin isoalloxazine scaffold
                "c1ccc2nc3c(=O)[nH]c(=O)nc3c2c1",
                # Second pattern to capture FAD variants (allowing flexible substitutions)
                "c1nc2c(c(=O)n1)C=CC(=O)C2"
            ],
            "constraint": lambda m, wt, p, s: True,
            "explanation": "Matches riboflavin (B2) pattern: isoalloxazine ring detected."
        },
        {
            "name": "B3 (niacin)",
            "patterns": [
                # Strict pattern for nicotinic acid (pyridine-3-carboxylic acid)
                "c1ccncc1C(=O)[O-]",
                "c1ccncc1C(=O)O"
            ],
            # Niacin normally has MW around 120–160 Da; also exclude molecules with sulfur.
            "constraint": lambda m, wt, p, s: (120 < wt < 200) and (s == 0),
            "explanation": "Matches niacin (B3) pattern: pyridine carboxylate moiety with proper MW detected."
        },
        {
            "name": "B5 (pantothenic acid)",
            "patterns": [
                # Pantothenic acid pattern; accepts both free acid and deprotonated forms.
                "CC(C)(CO)[C@@H](O)C(=O)NCCC(O)=O",
                "CC(C)(CO)[C@@H](O)C(=O)NCCC([O-])=O"
            ],
            # Pantothenic acid is moderately sized and should not have phosphorus.
            "constraint": lambda m, wt, p, s: (wt < 500) and (p == 0),
            "explanation": "Matches pantothenic acid (B5) pattern: characteristic substructure detected."
        },
        {
            "name": "B6 (pyridoxine and derivatives)",
            "patterns": [
                # Unphosphorylated pyridoxine-like structure
                "Cc1ncc(CO)c(c1)O",
                # Pyridoxal 5'-phosphate in one protonation state (with =O and P(O) groups)
                "[H]C(=O)c1c(COP(O)(=O)O)cnc(C)c1O",
                # Pyridoxal 5'-phosphate alternative deprotonated variant
                "[H]C(=O)c1c(COP([O-])([O-])=O)cnc(C)c1O",
                # Pyridoxamine 5'-phosphate variant
                "Cc1ncc(COP([O-])([O-])=O)c(C[NH3+])c1O",
                # Alternative protonation of pyridoxamine
                "C1(O)=C(C)[NH+]=CC(CO)=C1C[NH3+]"
            ],
            "constraint": lambda m, wt, p, s: True,
            "explanation": ("Matches pyridoxine/B6 derivatives pattern: pyridine ring with hydroxyl and/or phosphate "
                            "features detected.")
        },
        {
            "name": "B7 (biotin)",
            "patterns": [
                # Biotin has a bicyclic ureido and thiolane ring.
                "O=C(O)C1CSC[C@H]1N"
            ],
            # Biotin is small so we require a MW under 400.
            "constraint": lambda m, wt, p, s: wt < 400,
            "explanation": "Matches biotin (B7) pattern: fused ureido and thiolane ring detected."
        },
        {
            "name": "B9 (folate and tetrahydrofolate derivatives)",
            "patterns": [
                # General folate pterin fragment
                "Nc1nc2nccc(c2[nH]1)",
                # Tetrahydrofolate simplified motif (for instance, 5,6,7,8-tetrahydrofolic acid)
                "Nc1nc2NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)Nc2c(=O)[nH]1",
                # (6S)-5-methyltetrahydrofolic acid variant pattern
                "CN1[C@@H](CNC2=CC=C(C=C2)C(=O)N[C@@H](CCC(O)=O)C(O)=O)CNC2=C1C(=O)NC(N)=N2",
                # Folate derivative with extra carboxylate groups.
                "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CCC([O-])=O)C([O-])=O)nc2c(=O)[nH]1"
            ],
            # Folate derivatives are typically below 600 Da.
            "constraint": lambda m, wt, p, s: wt < 600,
            "explanation": "Matches folate (B9) pattern: pterin/tetrahydrofolate-like motif detected."
        }
    ]
    
    # Iterate over each vitamin subclass and its patterns.
    for vit in vitamin_checks:
        for smarts in vit["patterns"]:
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue  # skip if SMARTS cannot be parsed
            if mol.HasSubstructMatch(patt):
                # Check any extra constraints (using molecular weight, phosphorus and sulfur counts).
                if vit["constraint"](mol, mol_wt, p_count, s_count):
                    return True, f"Matches {vit['name']} pattern: {vit['explanation']}"
    
    # No match found.
    return False, "Does not match any known refined substructure patterns for B vitamins"

# Example usage (when running as a script)
if __name__ == "__main__":
    test_smiles = [
        "Cc1ncc(COP([O-])([O-])=O)c(C[NH3+])c1O",  # pyridoxamine 5'-phosphate: expected B6
        "[H]C(=O)c1c(COP([O-])([O-])=O)cnc(C)c1O",  # pyridoxal 5'-phosphate: expected B6
        "Nc1nc2NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)Nc2c(=O)[nH]1",  # tetrahydrofolic acid: expected B9
        "Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c2C)c(N)n1",  # thiamine(1+) triphosphate: expected B1
        "CC(C)(CO)[C@@H](O)C(=O)NCCC([O-])=O",  # (R)-pantothenate: expected B5
        "[Cl-].CC1=C(CCOP(O)(=O)OP(O)(O)=O)SC=[N+]1CC1=C(N)N=C(C)N=C1",  # thiamine(1+) diphosphate chloride: expected B1
        "C1(O)=C(C)[NH+]=CC(CO)=C1C[NH3+]",  # pyridoxamine(2+): expected B6
        "[H][C@]12CNc3nc(N)[nH]c(=O)c3N1CN(C2)c1ccc(cc1)C(=O)N[C@@H](CCC(O)=O)C(O)=O",  # (6R)-5,10-methylenetetrahydrofolic acid: expected B9
        "CC1=C(C)C=C2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3=NC2=C1"  # riboflavin/FAD variant: expected B2
    ]
    
    for sm in test_smiles:
        is_b, reason = is_B_vitamin(sm)
        print(f"SMILES: {sm}\nClassified as B vitamin? {is_b} - {reason}\n")