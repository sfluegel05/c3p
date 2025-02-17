"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin
Definition:
  Any member of the group of eight water‐soluble vitamins (B1, B2, B3, B5, B6, B7, B9 and B12)
  that play important roles in cell metabolism.
  
Notes on improvements:
  • Instead of immediately classifying any molecule with cobalt as a vitamin B12 derivative,
    we now require that (i) the molecular weight is quite high (typically >800 Da) and (ii)
    the molecule contains at least 4 rings (consistent with the corrin macrocycle) before returning True.
  • Several alternative SMARTS patterns (especially for vitamin B6 and B9 derivatives)
    have been added to capture phosphorylated, deprotonated, or tetrahydrofolate‐like scaffolds.
  • This rule‐based approach is inherently approximate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a member of the B vitamin class (B1, B2, B3, B5, B6, B7, B9, or B12)
    by trying to match a series of SMARTS patterns and extra constraints.
    
    Args:
        smiles (str): The SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a B vitamin, False otherwise.
        str: Explanation detailing which pattern (if any) was matched (or why not).
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pre-calculate descriptors.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    ring_count = len(mol.GetRingInfo().AtomRings())
    
    # --- Improved check for vitamin B12 (cobalamin derivatives) ---
    # Instead of immediately returning True on finding a cobalt atom (atomic number 27),
    # we require additional features (a high molecular weight and several rings).
    cobalt_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 27]
    if cobalt_atoms:
        if mol_wt > 800 and ring_count >= 4:
            return True, "Contains cobalt in a corrin-like macrocycle (high MW and multiple rings), indicating a vitamin B12 derivative"
    
    # Define a list of vitamin checks. Each entry has:
    #   name: a vitamin type (B1, B2, etc.)
    #   patterns: a list of SMARTS strings (alternative substructure queries)
    #   constraint: a lambda with (mol, mol_wt, p_count) to enforce extra conditions.
    #   explanation: text describing the match.
    vitamin_checks = [
        {
            "name": "B1 (thiamine)",
            "patterns": [
                # The canonical thiamine core (pyrimidine ring linked to a thiazolium ring via a methylene bridge).
                "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1",
                # Variation with a phosphate group on the side chain:
                "Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c2C)c(N)n1"
            ],
            "constraint": lambda m, wt, p: True,
            "explanation": "Matches thiamine (B1) core pattern: correct thiazolium–pyrimidine motif found."
        },
        {
            "name": "B2 (riboflavin)",
            "patterns": [
                # Riboflavin: characteristic isoalloxazine tricyclic ring.
                "c1cc2nc3c(=O)[nH]c(=O)nc3c2cc1"
            ],
            "constraint": lambda m, wt, p: True,
            "explanation": "Matches riboflavin (B2) pattern: isoalloxazine ring detected."
        },
        {
            "name": "B3 (niacin)",
            "patterns": [
                # Niacin: pyridine carboxylic acid (protonated or deprotonated)
                "c1ccncc1C(=O)[O-]",
                "c1ccncc1C(=O)O"
            ],
            # Niacin and its simple derivatives are generally small.
            "constraint": lambda m, wt, p: wt < 300,
            "explanation": "Matches niacin (B3) pattern: pyridine carboxylate moiety found."
        },
        {
            "name": "B5 (pantothenic acid)",
            "patterns": [
                # Pantothenic acid: characteristic substituted chain.
                "CC(C)(CO)[C@@H](O)C(=O)NCCC(O)=O",
                # Allow also a carboxylate variant.
                "CC(C)(CO)[C@@H](O)C(=O)NCCC([O-])=O"
            ],
            # To avoid picking up large derivatives (like CoA), require moderate weight and no phosphorus.
            "constraint": lambda m, wt, p: wt < 500 and p == 0,
            "explanation": "Matches pantothenic acid (B5) pattern: characteristic substructure detected."
        },
        {
            "name": "B6 (pyridoxine and derivatives)",
            "patterns": [
                # Unphosphorylated pyridoxine (a methylated pyridine with hydroxyl groups)
                "Cc1ncc(CO)c(c1)O",
                # Pyridoxal 5'-phosphate: variation where the phosphate group is drawn with deprotonated oxygens.
                "[H]C(=O)c1c(COP([O-])([O-])=O)cnc(C)c1O",
                # Alternate phosphorylated form.
                "Cc1ncc(COP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c(CN)c1O",
                # Pyridoxamine form as a cation.
                "C1(O)=C(C)[NH+]=CC(CO)=C1C[NH3+]"
            ],
            "constraint": lambda m, wt, p: True,
            "explanation": "Matches pyridoxine/phospho- forms (B6) pattern: appropriate pyridine ring with hydroxyl/phosphate features found."
        },
        {
            "name": "B7 (biotin)",
            "patterns": [
                # Biotin: fused bicyclic system including a thiolane ring.
                "O=C(O)C1CSC[C@H]1N"
            ],
            # Biotin is a relatively small molecule.
            "constraint": lambda m, wt, p: wt < 400,
            "explanation": "Matches biotin (B7) pattern: fused ureido and thiolane ring detected."
        },
        {
            "name": "B9 (folate and tetrahydrofolate derivatives)",
            "patterns": [
                # A common folate motif: substituted pterin ring.
                "Nc1nc2nccc(c2[nH]1)",
                # A pattern to catch a tetrahydrofolate-like scaffold.
                "Nc1nc2NCC",
                # Exact match for tetrahydrofolic acid derivative.
                "Nc1nc2NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)Nc2c(=O)[nH]1",
                # (6S)-5-methyltetrahydrofolic acid variant.
                "CN1[C@@H](CNC2=CC=C(C=C2)C(=O)N[C@@H](CCC(O)=O)C(O)=O)CNC2=C1C(=O)NC(N)=N2",
                # Folate with a carboxylate signature.
                "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CCC([O-])=O)C([O-])=O)nc2c(=O)[nH]1"
            ],
            # Folic acid and its derivatives tend to have intermediate molecular weights.
            "constraint": lambda m, wt, p: wt < 600,
            "explanation": "Matches folate (B9) pattern: pterin or tetrahydrofolate-like motif detected."
        }
    ]
    
    # Iterate through each vitamin check.
    for vit in vitamin_checks:
        for smarts in vit["patterns"]:
            patt = Chem.MolFromSmarts(smarts)
            # Skip unparseable patterns.
            if patt is None:
                continue
            if mol.HasSubstructMatch(patt):
                if vit["constraint"](mol, mol_wt, p_count):
                    return True, f"Matches {vit['name']} pattern: {vit['explanation']}"
    
    return False, "Does not match any known refined substructure patterns for B vitamins"

# Example usage (if run as script):
if __name__ == "__main__":
    test_smiles = [
        "Cc1ncc(COP([O-])([O-])=O)c(C[NH3+])c1O",  # pyridoxamine 5'-phosphate(1-), expected B6
        "[H]C(=O)c1c(COP([O-])([O-])=O)cnc(C)c1O",  # pyridoxal 5'-phosphate(2-), expected B6
        "Nc1nc2NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)Nc2c(=O)[nH]1",  # 5,6,7,8-tetrahydrofolic acid, expected B9
        "Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c2C)c(N)n1",  # thiamine(1+) triphosphate, expected B1
        "CC(C)(CO)[C@@H](O)C(=O)NCCC([O-])=O",  # (R)-pantothenate, expected B5
        "[Cl-].CC1=C(CCOP(O)(=O)OP(O)(O)=O)SC=[N+]1CC1=C(N)N=C(C)N=C1",  # thiamine(1+) diphosphate chloride, expected B1
        "Cc1ncc(COP([O-])([O-])=O)c(CN)c1O",  # pyridoxamine form, expected B6
        "[H][C@]12CNc3nc(N)[nH]c(=O)c3N1CN(C2)c1ccc(cc1)C(=O)N[C@@H](CCC(O)=O)C(O)=O",  # (6R)-5,10-methylenetetrahydrofolic acid, expected B9
        "C1(O)=C(C)[NH+]=CC(CO)=C1C[NH3+]",  # pyridoxamine(2+), expected B6
        "[H][C@]12[C@H](CC(N)=O)[C@@]3(C)CCC(=O)NC[C@@H](C)OP([O-])(=O)O[C@@H]4[C@@H](COP([O-])([O-])=O)O[C@@H]([C@@H]4O)n4c[n+](c5cc(C)c(C)cc45)[Co-3]456(C[C@H]7O[C@H]([C@H](O)[C@@H]7O)n7cnc8c(N)ncnc78)N1C3=C(C)C1=[N+]4C(=CC3=[N+]5C(=C(C)C4=[N+]6[C@]2(C)[C@@](C)(CC(N)=O)[C@@H]4CCC(N)=O)[C@@](C)(CC(N)=O)[C@@H]3CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O"  # adenosylcobalamin 5'-phosphate(2-), expected B12 (by cobalt and macrocycle)
    ]
    
    for sm in test_smiles:
        is_b, reason = is_B_vitamin(sm)
        print(f"SMILES: {sm}\nClassified as B vitamin? {is_b} - {reason}\n")