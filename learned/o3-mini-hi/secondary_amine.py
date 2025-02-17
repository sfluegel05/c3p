"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: secondary amine 
Definition: A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.

This function uses two SMARTS patterns:
  1. A standard secondary amine (R2NH): the nitrogen is trigonal (three substituents) with one hydrogen
     and is attached to exactly two carbon (hydrocarbyl) groups. We exclude those nitrogens that are bound
     to a carbonyl (i.e. are part of an amide).
  2. A nitrosated derivative in which one of the hydrogen positions is taken by a nitroso (–N=O) group.
     Here the nitrogen has no hydrogen (H0) but three substituents (two carbons and one nitroso substituent).
For either case the effective hydrogen count (H + (nitroso_count)) is 1.
"""

from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if the molecule (given by its SMILES string)
    contains at least one secondary amine group as defined below:
      • It is derived from ammonia by replacing precisely two hydrogens with hydrocarbyl groups.
      • Optionally one hydrogen may be “replaced” by a nitroso (–N=O) group – in which case that group is counted as a hydrogen.
      • The candidate nitrogen must not be attached to a carbonyl group (to avoid amide N).
    
    The function attempts two SMARTS substructure searches:
       pat_standard: matches typical secondary amines ([NX3;H1;!$(N[C]=O)]([#6])([#6]))
       pat_nitroso: matches a nitrosated secondary amine ([NX3;H0;!$(N[C]=O)]([#6])([#6])([NX2]=O))
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if at least one secondary amine is found, False otherwise.
       str: An explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for a standard secondary amine: 
    # nitrogen (N) that is trivalent (NX3) with one attached hydrogen (H1),
    # NOT attached to a carbonyl carbon (the negative pattern !$(N[C]=O))
    # and having two substituents that are carbon ([#6]).
    pat_standard = Chem.MolFromSmarts("[NX3;H1;!$(N[C]=O)]([#6])([#6])")

    # SMARTS for a nitrosated secondary amine:
    # This pattern accepts a nitrogen with no attached hydrogen (H0) but with three substituents,
    # two of which are carbon and one is a nitroso group. In nitroso,
    # the substituent is defined as a nitrogen (NX2) that is double-bonded to an oxygen.
    pat_nitroso = Chem.MolFromSmarts("[NX3;H0;!$(N[C]=O)]([#6])([#6])([NX2]=O)")

    if mol.HasSubstructMatch(pat_standard):
        return True, ("Found a secondary amine substructure: a nitrogen with one hydrogen and two "
                      "hydrocarbyl (carbon) substituents that is not directly bonded to a carbonyl group.")
    elif mol.HasSubstructMatch(pat_nitroso):
        return True, ("Found a nitrosated secondary amine substructure: a nitrogen with two hydrocarbyl "
                      "substituents and a nitroso group (counted as replacing a hydrogen) and no direct "
                      "bond to a carbonyl group.")
    else:
        return False, ("No secondary amine detected: did not find a nitrogen atom with exactly two hydrocarbyl "
                       "substituents and one effective hydrogen (or nitroso) that is not acylated.")

# For debugging and testing the function with several examples
if __name__ == "__main__":
    # List of test SMILES as provided (only a few examples; the full list is long)
    test_smiles = [
        "C=1(C(=CC=C(C1)CCN[C@@H](CCC=2C=CC(=CC2)O)C)O)O",  # (R)-dobutamine
        "CNC1CCCCC1",  # N-methylcyclohexylamine
        "CC(C)NCC(C)(C)N",  # N(1)-isopropyl-2-methylpropan-1,2-diamine
        "CNc1ccccc1",  # N-methylaniline
        "O=NN1CCCCC1",  # N-nitrosopiperidine
        "[H]N(C)C",  # dimethylamine
        # A peptide fragment (should not be classified as a secondary amine)
        "O=C(NCC(=O)N[C@@H](CO)C(O)=O)[C@@H](N)CC=1NC=NC1"
    ]
    for smi in test_smiles:
        result, reason = is_secondary_amine(smi)
        print(f"SMILES: {smi}\n  -> {result}, Reason: {reason}\n")