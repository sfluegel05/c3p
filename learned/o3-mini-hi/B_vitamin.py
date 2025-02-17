"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin (vitamins B1, B2, B3, B5, B6, B7, B9 and B12)
Definition:
   Any member of the group of eight water‐soluble vitamins originally thought to be a single compound (vitamin B)
   that play important roles in cell metabolism.
   (Members: vitamin B1, B2, B3, B5, B6, B7, B9, and B12.)
   
Note:
  This implementation uses a series of SMARTS substructure searches that attempt to capture key structural
  motifs of each vitamin type. Because B vitamins are a chemically diverse group, this approach is approximate.
"""

from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) is a member of the B vitamin class.
    The function uses substructure matching. It first checks if cobalt is present (indicating B12)
    and then looks for patterns corresponding approximately to vitamins B1, B2, B3, B5, B6, B7, and B9.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a B vitamin, False otherwise.
        str: A short explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cobalt atom (atomic number 27) which is characteristic of vitamin B12 derivatives
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 27:
            return True, "Contains cobalt atom, indicating a vitamin B12 derivative"

    # Define a list of (vitamin name, SMARTS pattern, explanation) for other types of B vitamins
    vitamin_patterns = [
        # Vitamin B1 (thiamine): contains a thiazole ring fused to a pyrimidine.
        ("B1 (thiamine)", "c1ncsc1", "Contains thiazole ring typical of thiamine (B1)"),
        # Vitamin B2 (riboflavin): has an isoalloxazine ring.
        ("B2 (riboflavin)", "c1cc2nc3c(=O)[nH]c(=O)nc3c2cc1", "Contains isoalloxazine motif typical of riboflavin (B2)"),
        # Vitamin B3 (niacin): key motif is a pyridine ring usually substituted by a carboxyl group.
        ("B3 (niacin)", "c1ccncc1C(=O)[O-]", "Contains pyridine carboxylic acid moiety typical of niacin (B3)"),
        # Vitamin B5 (pantothenic acid): contains a substituted 2,4-dihydroxy acid amide fragment.
        ("B5 (pantothenic acid)", "CC(C)(CO)[C@@H](O)C(=O)NCCC", "Matches substructure typical of pantothenic acid (B5)"),
        # Vitamin B6 (pyridoxine): typically possesses a methylated pyridine ring with hydroxyl substituents.
        ("B6 (pyridoxine)", "Cc1ncc(CO)c(c1)O", "Contains pyridoxine-like ring common to vitamin B6 compounds"),
        # Vitamin B7 (biotin): characterized by a fused bicyclic ureido and tetrahydrothiophane ring.
        ("B7 (biotin)", "O=C(N[C@H]1CSC[C@H]1)O", "Matches partial biotin substructure (B7)"),
        # Vitamin B9 (folate): contains a pterin ring often connected to a para-aminobenzoate moiety.
        ("B9 (folate)", "Nc1nc2nccc(c2[nH]1)", "Contains pterin ring system common to folate derivatives (B9)")
    ]

    # Try matching each of the patterns
    for vit_name, smarts, explanation in vitamin_patterns:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            # If pattern cannot be parsed, skip it.
            continue
        if mol.HasSubstructMatch(patt):
            return True, f"Matches {vit_name} pattern: {explanation}"

    # If none of the patterns match, then classify as not a B vitamin.
    return False, "Does not match any known substructure patterns for B vitamins"

# Example usage:
if __name__ == "__main__":
    # Test a few examples
    test_smiles = [
        "Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c2C)c(N)n1",  # thiamine triphosphate (B1 variant)
        "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C[C@H](O)[C@H](O)[C@H](O)CO)C2cc1C",  # riboflavin (B2)
        "OC(=O)C(=O)c1ccncc1",  # a pyridine dicarboxylic acid (B3–like)
        "CC(C)(CO)[C@@H](O)C(=O)NCCC(O)=O",  # pantothenic acid (B5)
        "Cc1ncc(CO)c(c1)O",  # pyridoxine (B6)
        "O=C(N[C@H]1CSC[C@H]1)O",  # biotin partial (B7)
        "Nc1nc2nccc(c2[nH]1)"  # folate fragment (B9)
    ]
    
    for sm in test_smiles:
        is_B, reason = is_B_vitamin(sm)
        print(f"SMILES: {sm}\nClassified as B vitamin? {is_B} - {reason}\n")