"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: secondary amine 
Definition: A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
Only exocyclic secondary amine groups are counted (i.e. the N must not be part of a ring).
Note: An optional nitroso substitution (–N=O) is allowed (counted as a replacement for hydrogen).

The function attempts three SMARTS-based substructure searches:
  1. Standard (aliphatic) secondary amine:
      - Non-aromatic nitrogen with three neighbors (one hydrogen, two carbons)
      - Not directly bonded to a carbonyl carbon.
  2. Standard aromatic secondary amine:
      - An aromatic nitrogen drawn as lower-case [nH] that is bonded to two carbon atoms.
      - We again require that the nitrogen is exocyclic (not in a ring).
  3. Nitrosated secondary amine:
      - Non-aromatic nitrogen that has no H but three substituents (two carbons and one nitroso group)
      - Not bonded to a carbonyl.
After matching, we filter out any nitrogen that is a member of any ring.
"""

from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if the molecule (given by its SMILES string)
    contains at least one exocyclic secondary amine group as defined below:
      • Derived from ammonia by replacing exactly two hydrogens with hydrocarbyl groups.
      • Optionally one hydrogen is “replaced” by a nitroso (–N=O) group (this counts as a hydrogen).
      • The candidate nitrogen must not be attached to a carbonyl group.
      • The nitrogen must not be part of a ring.
      
    The function uses three SMARTS patterns:
       pat_standard_aliphatic: typical secondary amine with one hydrogen (non‐aromatic)
       pat_standard_aromatic: an aromatic secondary amine (e.g. from pyrrole‐like depiction, though in our case we exclude ring members)
       pat_nitroso: secondary amine where one hydrogen is replaced by a nitroso group.

    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if an exocyclic secondary amine is found, False otherwise.
       str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for a standard secondary amine (aliphatic, non-aromatic):
    # [NX3;H1;!$(N[C]=O)] ensures a trigonal (sp3) nitrogen with one hydrogen and not bound to C=O.
    # ([#6])([#6]) adds requirement for two carbon substituents.
    pat_standard_aliphatic = Chem.MolFromSmarts("[NX3;H1;!$(N[C]=O)]([#6])([#6])")
    
    # SMARTS for a standard secondary amine (aromatic variant):
    # [nH] matches an aromatic nitrogen with one hydrogen.
    # We require two carbon neighbors.
    pat_standard_aromatic = Chem.MolFromSmarts("[nH]([#6])([#6])")

    # SMARTS for nitrosated secondary amine:
    # [NX3;H0;!$(N[C]=O)] matches a (likely non-aromatic) nitrogen with no hydrogen,
    # but with two carbons and one nitroso group ([NX2]=O) as substituents.
    pat_nitroso = Chem.MolFromSmarts("[NX3;H0;!$(N[C]=O)]([#6])([#6])([NX2]=O)")

    # Function to filter out matches where the matched nitrogen is part of a ring.
    def valid_match(match, smarts_type=''):
        # In our SMARTS, the first atom is the candidate nitrogen.
        n_idx = match[0]
        atom = mol.GetAtomWithIdx(n_idx)
        # Reject if the nitrogen is in any ring.
        if atom.IsInRing():
            return False
        return True

    # Search for aliphatic secondary amine matches.
    matches_aliphatic = mol.GetSubstructMatches(pat_standard_aliphatic)
    for match in matches_aliphatic:
        if valid_match(match, 'aliphatic'):
            return True, ("Found an exocyclic secondary amine substructure (aliphatic variant): "
                          "a nitrogen with one hydrogen and two hydrocarbyl (carbon) substituents that is not "
                          "directly bonded to a carbonyl group and not part of a ring.")

    # Search for aromatic secondary amine matches.
    matches_aromatic = mol.GetSubstructMatches(pat_standard_aromatic)
    for match in matches_aromatic:
        if valid_match(match, 'aromatic'):
            return True, ("Found an exocyclic secondary amine substructure (aromatic variant): "
                          "an aromatic nitrogen with one hydrogen and two hydrocarbyl (carbon) substituents, not "
                          "bonded to a carbonyl group and not within a ring system.")

    # Search for nitrosated secondary amine matches.
    matches_nitroso = mol.GetSubstructMatches(pat_nitroso)
    for match in matches_nitroso:
        if valid_match(match, 'nitroso'):
            return True, ("Found a nitrosated secondary amine substructure: a nitrogen with two hydrocarbyl "
                          "substituents and a nitroso group (counted as replacing a hydrogen) that is not directly "
                          "bonded to a carbonyl group and is exocyclic.")

    return False, ("No exocyclic secondary amine detected: did not find a nitrogen atom with exactly two hydrocarbyl "
                   "substituents and one effective hydrogen (or replaced by a nitroso group), or all candidate nitrogen(s) "
                   "were part of ring systems.")

# For debugging and testing the function with a few examples:
if __name__ == "__main__":
    # A selected list of test SMILES (not the whole list provided)
    test_molecules = [
        ("(R)-dobutamine", "C=1(C(=CC=C(C1)CCN[C@@H](CCC=2C=CC(=CC2)O)C)O)O"),
        ("N-methylcyclohexylamine", "CNC1CCCCC1"),
        ("N(1)-isopropyl-2-methylpropan-1,2-diamine", "CC(C)NCC(C)(C)N"),
        ("(S)-dobutamine", "C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1"),
        ("N-(3-aminopropyl)-4-aminobutanal", "[H]C(=O)CCCNCCCN"),
        ("bisoprolol", "CC(C)NCC(O)COc1ccc(COCCOC(C)C)cc1"),
        ("N-methylaniline", "CNc1ccccc1"),
        ("nortriptyline", "CNCCC=C1c2ccccc2CCc2ccccc12"),
        ("1H-pyrrole", "c1cc[nH]c1"),  # This one should NOT be classified (nitrogen in ring)
        ("pyrrolidine", "C1CCNC1"),      # This one should NOT be classified (nitrogen in ring)
        ("dimethylamine", "[H]N(C)C"),
    ]
    
    for name, smi in test_molecules:
        result, explanation = is_secondary_amine(smi)
        print(f"SMILES: {smi}\n  NAME: {name}\n  -> {result}, Reason: {explanation}\n")