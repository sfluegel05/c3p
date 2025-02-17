"""
Classifies: CHEBI:17478 aldehyde
"""
#!/usr/bin/env python
"""
Classifies: Aldehyde – a compound having the functional group RC(=O)H 
(including formaldehyde H2C=O) where the carbonyl carbon is bonded 
to exactly one hydrogen and one R group (or two hydrogens for formaldehyde).

This version improves on a previous attempt by explicitly using SMARTS patterns
that capture typical aldehyde connectivity. For a typical (non‐formaldehyde)
aldehyde we require that the carbonyl carbon (which carries one hydrogen)
is attached to two carbon atoms (one coming from the R group and the other
as part of the aldehyde group binding) and that it is exocyclic. For formaldehyde,
the pattern is adapted to require two hydrogens.

False positives (often from lactones or esters) are reduced by checking that
the candidate carbonyl carbon (from the SMARTS match) is not part of any ring.
"""

from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    
    An aldehyde is defined as a compound containing a carbonyl group (C=O)
    in which the carbonyl carbon is bonded to exactly one hydrogen and one R group
    (or two hydrogens in the case of formaldehyde). In addition we require that
    the carbonyl carbon be exocyclic (i.e. not in a ring) to avoid matching cyclic 
    acyl groups (e.g. in lactones) that have similar connectivity.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if at least one valid aldehyde group is found, False otherwise.
        str: Explanation for classification.
    """
    # Parse SMILES and add explicit hydrogens so connectivity is unambiguous.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define two SMARTS patterns:
    # 1. For a typical aldehyde: the pattern "[#6][C;H1](=O)[#6]" means:
    #    - a carbon (atomic number 6) attached to a carbonyl carbon that has exactly one hydrogen (H1)
    #      and double bonded to an oxygen, followed by a carbon substituent.
    #
    # 2. For formaldehyde: the pattern "[C;H2](=O)" means:
    #    - a carbon (without an R group) carrying two hydrogens, double-bonded to an oxygen.
    patt_aldehyde = Chem.MolFromSmarts("[#6][C;H1](=O)[#6]")
    patt_formaldehyde = Chem.MolFromSmarts("[C;H2](=O)")
    
    if patt_aldehyde is None or patt_formaldehyde is None:
        return False, "Failed to compile SMARTS patterns"
    
    valid_matches = 0
    
    # Check typical aldehyde pattern matches.
    matches = mol.GetSubstructMatches(patt_aldehyde)
    for match in matches:
        # In the pattern "[#6][C;H1](=O)[#6]" the second atom is the carbonyl carbon.
        carbonyl_idx = match[1]
        carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
        # Only accept if the carbonyl carbon is exocyclic (not in a ring).
        if carbonyl.IsInRing():
            continue
        valid_matches += 1

    # Check formaldehyde pattern matches.
    matches = mol.GetSubstructMatches(patt_formaldehyde)
    for match in matches:
        # In the pattern "[C;H2](=O)" the first atom is the carbonyl carbon.
        carbonyl_idx = match[0]
        carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
        # Accept only if this carbon is not in a ring.
        if carbonyl.IsInRing():
            continue
        valid_matches += 1

    if valid_matches == 0:
        return False, "No valid aldehyde group found (check connectivity and ensure exocyclic C=O with required H count)"
    
    # Form the explanation
    group_word = "group" if valid_matches == 1 else "groups"
    reason = f"Found {valid_matches} aldehyde {group_word} with proper connectivity (exocyclic C=O with one H and one R group, or two H's for formaldehyde)"
    return True, reason

# Optional testing code
if __name__ == "__main__":
    # List of some example SMILES strings provided in the prompt.
    test_smiles = [
        "O=CC(CCC=C(C)C)C",  # 5-Heptenal, 2,6-dimethyl-
        "Oc1c(C=O)ccc2ccccc12",  # 1-hydroxy-2-naphthaldehyde
        "CC(=O)O[C@H]1CC[C@]2(C=O)[C@H]3CC[C@@]4(C)[C@@H](CCC4=O)[C@@H]3CC=C2C1",  # 3beta-Hydroxy-17-oxoandrost-5-en-19-al acetate
        "C(CCCCCCCCCCC(=O)[H])CCCCC",  # heptadecanal
        "O=CCCCCCCCCCCCCCCCCCCCCCCC",  # tetracosanal
        "CCCC/C=C/C=O",  # (E)-hept-2-enal
        "CCCCCCCCCCCCCCCCC=O",  # octadecanal
        "[H]C(=O)c1cc(N)ccc1O",  # 5-aminosalicylaldehyde
        "[H]C(=O)c1cccc(O)c1O",  # 2,3-dihydroxybenzaldehyde
        "Cc1cc(O)c(C=O)c2Oc3cc4OC(C)(C)CC(=O)c4c(C)c3OC(=O)c12",  # mollicellin B
        "C[C@H](CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C)C=O",  # (25R)-3-ketocholest-4-en-26-al
        "OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC=O",  # 20-oxoarachidonic acid
        "C(C#C)=O",  # prop-2-ynal
        "[H]C(=O)",  # formaldehyde
        "CC(C)CC=O",  # 3-methylbutanal
        # Some examples that previously led to false positives:
        "O=C1OC2(O[C@H]3[C@@H]([C@]4(CCC5=C([C@@]4(C3)C)CCC6[C@@]5(CC[C@H]([C@]6(C=O)C)O)C)C)[C@@H](C2)C)C(=C1C)C",  # Hexagonin D (should not classify)
        "[H]C(=O)CCCC([O-])=O",  # 5-oxopentanoate (should not classify)
    ]
    
    for sm in test_smiles:
        result, explanation = is_aldehyde(sm)
        print(f"SMILES: {sm}\nResult: {result}\nExplanation: {explanation}\n")