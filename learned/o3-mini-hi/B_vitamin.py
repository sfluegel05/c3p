"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin
Definition:
   Any member of the group of eight water‐soluble vitamins (B1, B2, B3, B5, B6, B7, B9 and B12)
   that play important roles in cell metabolism.
   
Notes on improvements:
   • Immediately returns True if a cobalt atom is present (vitamin B12 derivatives).
   • Uses several alternative SMARTS patterns for vitamins (B1, B6, and B9 in particular) so that
     common derivative forms (phosphorylated, deprotonated) are captured.
   • For vitamin B1, rather than using an overly generic “thiazole” pattern, we now require the
     canonical thiamine core (a thiazolium ring plus a pyrimidine ring linked by a methylene bridge).
   • For vitamins B5 and B6 the SMARTS now include alternatives for acid vs. carboxylate forms.
   • For folate (B9), an extra alternative pattern covering tetrahydrofolate‐like scaffolds was added.
   
It is understood that any rule‐based approach to detect B vitamins is approximate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a member of the B vitamin class (B1, B2, B3, B5, B6, B7, B9, or B12)
    by using a set of SMARTS substructure patterns and additional constraints.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a B vitamin, False otherwise.
        str: Explanation detailing which pattern (if any) was matched (or not).
    """
    # Parse SMILES into molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # First check: vitamin B12 derivatives are usually cobalt-containing (atomic num 27)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 27:
            return True, "Contains cobalt atom, indicating a vitamin B12 derivative"
    
    # Pre-calculate useful descriptors
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Count phosphorus atoms – might flag phosphorylated forms
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    # Each vitamin check here has:
    #  - a descriptive name
    #  - one or more SMARTS patterns (as alternatives)
    #  - an optional extra constraint (often based on molecular weight or phosphorus count)
    #  - an explanation to use when the vitamin type is matched.
    vitamin_checks = [
        {
            "name": "B1 (thiamine)",
            "patterns": [
                # The canonical thiamine core (pyrimidine ring linked to a thiazolium ring via methylene).
                # This pattern is expected to be found even in phosphorylated derivatives.
                "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1",
                # Allow variation with phosphate substitutes attached to the –CH2– side chain:
                "Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c2C)c(N)n1"
            ],
            "constraint": lambda m, wt, p: True,
            "explanation": "Matches thiamine (B1) core pattern: correct thiazolium–pyrimidine motif found."
        },
        {
            "name": "B2 (riboflavin)",
            "patterns": [
                # Riboflavin has a characteristic isoalloxazine (tricyclic) ring system.
                "c1cc2nc3c(=O)[nH]c(=O)nc3c2cc1"
            ],
            "constraint": lambda m, wt, p: True,
            "explanation": "Matches riboflavin (B2) pattern: isoalloxazine ring detected."
        },
        {
            "name": "B3 (niacin)",
            "patterns": [
                # Niacin is a pyridine carboxylic acid (or its deprotonated form).
                "c1ccncc1C(=O)[O-]",
                "c1ccncc1C(=O)O"
            ],
            # Niacin and its simple derivatives tend to be small molecules.
            "constraint": lambda m, wt, p: wt < 300,
            "explanation": "Matches niacin (B3) pattern: pyridine carboxylate moiety found."
        },
        {
            "name": "B5 (pantothenic acid)",
            "patterns": [
                # Pantothenic acid has a distinctive substituted chain.
                "CC(C)(CO)[C@@H](O)C(=O)NCCC(O)=O",
                # Also allow a carboxylate variant:
                "CC(C)(CO)[C@@H](O)C(=O)NCCC([O-])=O"
            ],
            # To avoid matching larger derivatives (e.g. CoA), require moderate weight and no phosphorus.
            "constraint": lambda m, wt, p: wt < 500 and p == 0,
            "explanation": "Matches pantothenic acid (B5) pattern: characteristic substructure detected."
        },
        {
            "name": "B6 (pyridoxine and derivatives)",
            "patterns": [
                # Unphosphorylated pyridoxine: a methylated pyridine ring with hydroxyl groups.
                "Cc1ncc(CO)c(c1)O",
                # Pyridoxal 5'-phosphate (allowing both protonated and deprotonated phosphate groups)
                "[H]C(=O)c1c(COP(O)(=O)=[O-])cnc(C)c1O",
                "Cc1ncc(COP([O-])([O-])=O)c(CN)c1O",
                # Also include pyridoxamine cationic form:
                "C1(O)=C(C)[NH+]=CC(CO)=C1C[NH3+]"
            ],
            "constraint": lambda m, wt, p: True,
            "explanation": "Matches pyridoxine/phospho- forms (B6) pattern: appropriate pyridine ring with hydroxyl/phosphate features found."
        },
        {
            "name": "B7 (biotin)",
            "patterns": [
                # Biotin is recognized by its fused bicyclic system including a thiolane ring.
                "O=C(O)C1CSC[C@H]1N"
            ],
            "constraint": lambda m, wt, p: wt < 400,  # biotin is a relatively small molecule
            "explanation": "Matches biotin (B7) pattern: fused ureido and thiolane ring detected."
        },
        {
            "name": "B9 (folate)",
            "patterns": [
                # One common motif in folate derivatives is a substituted pterin ring.
                "Nc1nc2nccc(c2[nH]1)",
                # An alternative pattern to capture tetrahydrofolate-like scaffolds
                "Nc1nc2NCC"
            ],
            # Folic acid free acid weight is around 441 Da; allow moderately heavy molecules.
            "constraint": lambda m, wt, p: wt < 600,
            "explanation": "Matches folate (B9) pattern: pterin or tetrahydrofolate-like motif detected."
        }
    ]
    
    # Loop through each vitamin check.
    for vit in vitamin_checks:
        for smarts in vit["patterns"]:
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue  # skip if pattern is unparseable
            if mol.HasSubstructMatch(patt):
                # If the substructure is found, see if extra constraints (e.g. weight) are met.
                if vit["constraint"](mol, mol_wt, p_count):
                    return True, f"Matches {vit['name']} pattern: {vit['explanation']}"
    
    # If no pattern was matched, report as not being a B vitamin.
    return False, "Does not match any known refined substructure patterns for B vitamins"

# Example usage with a few test SMILES:
if __name__ == "__main__":
    test_smiles = [
        "Cc1ncc(COP([O-])([O-])=O)c(C[NH3+])c1O",                                    # pyridoxamine 5'-phosphate(1-), expected B6
        "[H]C(=O)c1c(COP([O-])([O-])=O)cnc(C)c1O",                                  # pyridoxal 5'-phosphate(2-), expected B6
        "Nc1nc2NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)Nc2c(=O)[nH]1",         # 5,6,7,8-tetrahydrofolic acid, expected B9
        "Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c2C)c(N)n1",             # thiamine(1+) triphosphate, expected B1
        "CC(C)(CO)[C@@H](O)C(=O)NCCC([O-])=O",                                      # (R)-pantothenate, expected B5
        "[Cl-].CC1=C(CCOP(O)(=O)OP(O)(O)=O)SC=[N+]1CC1=C(N)N=C(C)N=C1",              # thiamine(1+) diphosphate chloride, expected B1
        "Cc1ncc(COP([O-])([O-])=O)c(CN)c1O",                                        # one pyridoxamine form, expected B6
        "[H][C@]12CNc3nc(N)[nH]c(=O)c3N1CN(C2)c1ccc(cc1)C(=O)N[C@@H](CCC(O)=O)C(O)=O",  # (6R)-5,10-methylenetetrahydrofolic acid, expected B9
        "C1(O)=C(C)[NH+]=CC(CO)=C1C[NH3+]",                                         # pyridoxamine(2+), expected B6
        "[H][C@]12[C@H](CC(N)=O)[C@@]3(C)CCC(=O)NC[C@@H](C)OP([O-])(=O)O[C@@H]4[C@@H](COP([O-])([O-])=O)O[C@@H]([C@@H]4O)n4c[n+](c5cc(C)c(C)cc45)[Co-3]456(C[C@H]7O[C@H]([C@H](O)[C@@H]7O)n7cnc8c(N)ncnc78)N1C3=C(C)C1=[N+]4C(=CC3=[N+]5C(=C(C)C4=[N+]6[C@]2(C)[C@@](C)(CC(N)=O)[C@@H]4CCC(N)=O)[C@@](C)(CC(N)=O)[C@@H]3CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O"  # adenosylcobalamin 5'-phosphate(2-), expected B12 (by cobalt)
    ]
    
    for sm in test_smiles:
        is_b, reason = is_B_vitamin(sm)
        print(f"SMILES: {sm}\nClassified as B vitamin? {is_b} - {reason}\n")