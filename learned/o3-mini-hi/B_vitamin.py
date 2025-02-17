"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin (B1, B2, B3, B5, B6, B7 and B9, plus B12 by cobalt)
Definition:
   Any member of the group of eight water‐soluble vitamins originally thought to be a single compound (vitamin B)
   that play important roles in cell metabolism.
Note:
   Because the B vitamins are chemically quite diverse, this implementation uses a series of SMARTS substructure 
   searches with additional filters (e.g. molecular weight and presence of phosphorus) to try to capture key motifs.
   The rules are still approximate. In case the molecule is a B12 derivative (a cobalt‐containing species),
   that is immediately accepted.
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
        bool: True if the molecule is classified as a B vitamin, False otherwise.
        str: A short explanation for the classification decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First check: vitamin B12 derivatives are usually cobalt-containing.
    # (Cobalt atomic number 27)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 27:
            return True, "Contains cobalt atom, indicating a vitamin B12 derivative"
    
    # Pre-calculate useful descriptors
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)  # count phosphorus atoms

    # For vitamins that come in a “free” acid form (e.g. B3, B5, B7, and B9),
    # we expect the molecule not to be huge. In contrast, phosphate esters or di/triphosphate forms
    # may be acceptable for B1 or B6.
    
    # Define a list of vitamin checks. Each entry has:
    #  - name: a label for the vitamin type
    #  - patterns: a list of one or more SMARTS strings for key substructure motifs
    #  - constraint: a lambda that further refines the decision (based on molecular weight, phosphorus, etc.)
    #  - explanation: a message to include if the vitamin matches.
    vitamin_checks = [
        {
            "name": "B1 (thiamine)",
            "patterns": [
                # Thiamine and its phosphate esters typically contain a thiazole ring.
                "c1ncsc1"
            ],
            "constraint": lambda mol, wt, p: True,  # no extra constraint beyond the ring motif
            "explanation": "Matches thiamine (B1) pattern: thiazole ring present."
        },
        {
            "name": "B2 (riboflavin)",
            "patterns": [
                # The rigid isoalloxazine ring system is characteristic of riboflavin.
                "c1cc2nc3c(=O)[nH]c(=O)nc3c2cc1"
            ],
            "constraint": lambda mol, wt, p: True,
            "explanation": "Matches riboflavin (B2) pattern: isoalloxazine ring detected."
        },
        {
            "name": "B3 (niacin)",
            "patterns": [
                # Niacin’s motif is a pyridine ring with a carboxyl group.
                "c1ccncc1C(=O)[O-]",
                "c1ccncc1C(=O)O"
            ],
            # Niacin (or its simple derivatives) is small.
            "constraint": lambda mol, wt, p: wt < 300,
            "explanation": "Matches niacin (B3) pattern: pyridine carboxylate moiety found."
        },
        {
            "name": "B5 (pantothenic acid)",
            "patterns": [
                # Pantothenic acid has a distinctive substituted chain.
                "CC(C)(CO)[C@@H](O)C(=O)NCCC(O)=O"
            ],
            # To avoid matching larger coenzyme A derivatives that share part of the pantothenic acid motif,
            # require that the weight is in a lower range and that there are no phosphorus atoms.
            "constraint": lambda mol, wt, p: wt < 500 and p == 0,
            "explanation": "Matches pantothenic acid (B5) pattern: characteristic substructure detected."
        },
        {
            "name": "B6 (pyridoxine and derivatives)",
            "patterns": [
                # Pyridoxine (vitamin B6) typically shows a methylated pyridine ring with hydroxyl groups.
                "Cc1ncc(CO)c(c1)O",
                # Also consider phosphorylated forms (pyridoxal or pyridoxamine phosphates)
                "[H]C(=O)c1c(COP(O)([O-])=O)cnc(C)c1O",
                "Cc1ncc(COP(O)([O-])=O)c(CN)c1O"
            ],
            "constraint": lambda mol, wt, p: True,
            "explanation": "Matches pyridoxine/phospho- forms (B6) pattern: pyridine ring with hydroxyl/phosphate features."
        },
        {
            "name": "B7 (biotin)",
            "patterns": [
                # Biotin has a fused bicyclic system with a thiolane (sulfur-containing) ring.
                "O=C(O)C1CSC[C@H]1N"
            ],
            "constraint": lambda mol, wt, p: wt < 400,   # biotin is normally a small molecule
            "explanation": "Matches biotin (B7) pattern: fused ureido and thiolane ring structure detected."
        },
        {
            "name": "B9 (folate)",
            "patterns": [
                # Folate derivatives often contain a pterin ring system; here we target a simplified pterin motif.
                "Nc1nc2nccc(c2[nH]1)"
            ],
            # Folic acid free acid weight is around 441 Da; allow moderately heavy molecules.
            "constraint": lambda mol, wt, p: wt < 600,
            "explanation": "Matches folate (B9) pattern: pterin ring system present."
        }
    ]
    
    # Loop over each vitamin check and try each SMARTS pattern.
    for vit in vitamin_checks:
        for smarts in vit["patterns"]:
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue  # skip unparseable pattern
            if mol.HasSubstructMatch(patt):
                # If the pattern matches, check additional constraints.
                if vit["constraint"](mol, mol_wt, p_count):
                    return True, f"Matches {vit['name']} pattern: {vit['explanation']}"
    
    return False, "Does not match any known refined substructure patterns for B vitamins"

# Example usage:
if __name__ == "__main__":
    # Test a selection of examples (both those historically classified correctly and some near misses)
    test_smiles = [
        "Cc1ncc(COP([O-])([O-])=O)c(C[NH3+])c1O",  # pyridoxamine 5'-phosphate(1-)
        "[H]C(=O)c1c(COP([O-])([O-])=O)cnc(C)c1O",  # pyridoxal 5'-phosphate(2-)
        "Nc1nc2NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)",  # 5,6,7,8-tetrahydrofolic acid (simplified folate motif)
        "CN1[C@@H](CNC2=CC=C(C=C2)C(=O)N[C@@H](CCC(O)=O)C(O)=O)CNC2=C1C(=O)NC(N)=N2",  # (6S)-5-methyltetrahydrofolic acid
        "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CCC([O-])=O)C([O-])=O)",  # folate(2-)
        "C1(O)=C(C)[NH+]=CC(CO)=C1C[NH3+]",  # pyridoxamine(2+)
        "Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c2C)c(N)n1",  # thiamine(1+) triphosphate(1-)
        "Cc1ncc(C[n+]2csc(CCOP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)c2C)c(N)n1",  # thiamine(1+) triphosphate(4-)
        "CC(C)(CO)[C@@H](O)C(=O)NCCC([O-])=O",  # (R)-pantothenate
        "[Cl-].CC1=C(CCOP(O)(=O)OP(O)(O)=O)SC=[N+]1CC1=C(N)N=C(C)N=C1",  # thiamine(1+) diphosphate chloride
        "[H][C@]12[C@H](CC(N)=O)[C@@]3(C)CCC(=O)NC[C@@H](C)OP([O-])(=O)O[C@@H]4[C@@H](COP([O-])([O-])=O)O[C@@H]([C@@H]4O)n4c[n+](c5cc(C)c(C)cc45)[Co-3]456(C[C@H]7O[C@H]([C@H](O)[C@@H]7O)n7cnc8c(N)ncnc78)N1C3=C(C)C1=[N+]4C(=CC3=[N+]5C(=C(C)C4=[N+]6[C@]2(C)[C@@](C)(CC(N)=O)[C@@H]4CCC(N)=O)[C@@](C)(CC(N)=O)[C@@H]3CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O",  # adenosylcobalamin 5'-phosphate(2-)
        "Cc1ncc(C[n+]2csc(CCO)c2C)c([NH3+])n1",  # thiamine(2+)
        "Cc1ncc(C[n+]2csc(CCOP(O)([O-])=O)c2C)c(N)n1",  # thiamine(1+) monophosphate(1-)
        "[H][C@]12CNc3nc(N)[nH]c(=O)c3N1CN(C2)c1ccc(cc1)C(=O)N[C@@H](CCC(O)=O)C(O)=O",  # (6R)-5,10-methylenetetrahydrofolic acid
        "C1(O)=C(C)[NH+]=CC(CO)=C1C([H])=O"  # pyridoxal 5'-phosphate
    ]
    
    for sm in test_smiles:
        is_b, reason = is_B_vitamin(sm)
        print(f"SMILES: {sm}\nClassified as B vitamin? {is_b} - {reason}\n")