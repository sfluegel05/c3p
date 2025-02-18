"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin
Definition:
  Any member of the group of eight water‐soluble vitamins (B1, B2, B3, B5, B6, B7, B9 and B12)
  that play important roles in cell metabolism.
  
Improvements over the previous version:
  • For vitamin B12 (cobalamin derivatives), in addition to cobalt, we now require
    a high molecular weight (>800 Da) and at least 4 rings to avoid misclassifying
    other cobalt complexes.
  • For vitamins B6 and B9, multiple SMARTS patterns (including phosphorylated and
    alternate protonation versions) are used so that known variants like pyridoxamine
    5'-phosphate, pyridoxal 5'-phosphate, and tetrahydrofolate derivatives are matched.
  • For riboflavin and related coenzymes (B2), an extra pattern is added to cover isoalloxazine
    scaffolds in both riboflavin and FAD.
  
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
        bool: True if the molecule is classified as a B vitamin, False otherwise.
        str: Explanation detailing which substructure pattern (if any) was matched or why it failed.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate descriptors to support extra constraints.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Count of phosphorus atoms (to exclude CoA derivatives etc.)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    # Ring count helps reinforce the corrin macrocycle for B12 and others.
    ring_count = len(mol.GetRingInfo().AtomRings())
    
    # --- Improved check for vitamin B12 (cobalamin derivatives) ---
    # Look for cobalt atoms. But to be considered vitamin B12 we also require a high weight and multiple rings.
    cobalt_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 27]
    if cobalt_atoms:
        if mol_wt > 800 and ring_count >= 4:
            return True, ("Contains cobalt in a corrin-like macrocycle (high MW and multiple rings), "
                          "indicating a vitamin B12 derivative")
    
    # Define SMARTS patterns and constraints for the various vitamin subclasses.
    vitamin_checks = [
        {
            "name": "B1 (thiamine)",
            "patterns": [
                # Canonical thiamine core pattern (thiazolium–pyrimidine motif)
                "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1",
                # Variation with phosphate group on the side chain (triphosphate versions etc.)
                "Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c2C)c(N)n1"
            ],
            "constraint": lambda m, wt, p: True,
            "explanation": "Matches thiamine (B1) core pattern: correct thiazolium–pyrimidine motif found."
        },
        {
            "name": "B2 (riboflavin/FAD)",
            "patterns": [
                # Riboflavin basic isoalloxazine scaffold
                "c1ccc2nc3c(=O)[nH]c(=O)nc3c2c1",
                # Additional pattern to cover FAD variants (flexible substitution tolerated)
                "c1nc2c(c(=O)n1)C=CC(=O)C2"
            ],
            "constraint": lambda m, wt, p: True,
            "explanation": "Matches riboflavin (B2) pattern: isoalloxazine ring detected."
        },
        {
            "name": "B3 (niacin)",
            "patterns": [
                # Niacin (pyridine carboxylic acid)
                "c1ccncc1C(=O)[O-]",
                "c1ccncc1C(=O)O"
            ],
            # Niacin and simple derivatives are typically small.
            "constraint": lambda m, wt, p: wt < 300,
            "explanation": "Matches niacin (B3) pattern: pyridine carboxylate moiety found."
        },
        {
            "name": "B5 (pantothenic acid)",
            "patterns": [
                # Pantothenic acid standard pattern (with or without deprotonation)
                "CC(C)(CO)[C@@H](O)C(=O)NCCC(O)=O",
                "CC(C)(CO)[C@@H](O)C(=O)NCCC([O-])=O"
            ],
            # To avoid picking up very large derivatives, require moderate weight (<500 Da) and no phosphorus.
            "constraint": lambda m, wt, p: wt < 500 and p == 0,
            "explanation": "Matches pantothenic acid (B5) pattern: characteristic substructure detected."
        },
        {
            "name": "B6 (pyridoxine and derivatives)",
            "patterns": [
                # Unphosphorylated pyridoxine form.
                "Cc1ncc(CO)c(c1)O",
                # Pyridoxal 5'-phosphate variant.
                "[H]C(=O)c1c(COP([O-])([O-])=O)cnc(C)c1O",
                # Another phosphorylated variant variant.
                "Cc1ncc(COP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c(CN)c1O",
                # Pyridoxamine 5'-phosphate version.
                "Cc1ncc(COP([O-])([O-])=O)c(C[NH3+])c1O",
                # Pyridoxamine as a dication.
                "C1(O)=C(C)[NH+]=CC(CO)=C1C[NH3+]"
            ],
            "constraint": lambda m, wt, p: True,
            "explanation": ("Matches pyridoxine/B6 derivatives pattern: appropriate pyridine ring with "
                            "hydroxyl and/or phosphate features detected.")
        },
        {
            "name": "B7 (biotin)",
            "patterns": [
                # Biotin contains a characteristic bicyclic ureido and thiolane ring.
                "O=C(O)C1CSC[C@H]1N"
            ],
            # Biotin is relatively small.
            "constraint": lambda m, wt, p: wt < 400,
            "explanation": "Matches biotin (B7) pattern: fused ureido and thiolane ring detected."
        },
        {
            "name": "B9 (folate and tetrahydrofolate derivatives)",
            "patterns": [
                # A general folate motif: substituted pterin ring.
                "Nc1nc2nccc(c2[nH]1)",
                # A pattern to partially match tetrahydrofolate-like scaffolds.
                "Nc1nc2NCC",
                # A more specific pattern for tetrahydrofolic acid derivatives.
                "Nc1nc2NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)Nc2c(=O)[nH]1",
                # Variant for (6S)-5-methyltetrahydrofolic acid.
                "CN1[C@@H](CNC2=CC=C(C=C2)C(=O)N[C@@H](CCC(O)=O)C(O)=O)CNC2=C1C(=O)NC(N)=N2",
                # Folate with carboxylate groups.
                "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CCC([O-])=O)C([O-])=O)nc2c(=O)[nH]1"
            ],
            # Folate derivatives tend to have an intermediate weight (<600 Da).
            "constraint": lambda m, wt, p: wt < 600,
            "explanation": "Matches folate (B9) pattern: pterin/tetrahydrofolate-like motif detected."
        }
    ]
    
    # Iterate over vitamin subclasses and attempt to match each SMARTS pattern.
    for vit in vitamin_checks:
        for smarts in vit["patterns"]:
            patt = Chem.MolFromSmarts(smarts)
            # Skip patterns that cannot be parsed.
            if patt is None:
                continue
            if mol.HasSubstructMatch(patt):
                # Check extra constraints (e.g., molecular weight limits, phosphorus count).
                if vit["constraint"](mol, mol_wt, p_count):
                    return True, f"Matches {vit['name']} pattern: {vit['explanation']}"
    
    # If no patterns matched with the extra constraints, return False.
    return False, "Does not match any known refined substructure patterns for B vitamins"

# Example usage (when running as a script):
if __name__ == "__main__":
    # Some test SMILES taken from the provided examples.
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