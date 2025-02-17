"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: Phenylpropanoid – any aromatic compound based on a phenylpropane (C6–C3) skeleton.

Revised approach:
  • We define several SMARTS patterns capturing well‐known phenylpropanoid substructures:
      – Cinnamoyl ester/acid (as found in cinnamate esters/acids, isoeugenol, etc.)
      – Coumarin (lactone-fused benzene)
      – Flavone (a common flavonoid core)
      – Isoflavonoid (the isomer)
      – A generic unsaturated C6–C3 pattern: a benzene ring attached to an allylic (C=C) chain.
  • We removed the overly generic saturated C6–C3 pattern which was scoring many false positives.
  • If none of the above match but no benzene is even found, we classify as non-phenylpropanoid.
This heuristic still may not capture every design and may require further tuning.
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    
    Uses a set of SMARTS patterns that capture substructures typical of phenylpropanoid subclasses:
      • Cinnamoyl (acid or ester) moieties,
      • Coumarins (benzopyran-2-ones),
      • Flavone and isoflavonoid cores,
      • And a generic unsaturated C6–C3 chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phenylpropanoid, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns along with descriptive names.
    # The cinnamoyl patterns require a benzene ring attached via a trans double bond to a carbonyl (acid or ester)
    patterns = [
        ("Cinnamoyl Ester", "c1ccc(/C=C/C(=O)O)cc1"),
        ("Cinnamoyl Acid", "O=C(O)/C=C/c1ccccc1"),
        # Coumarin pattern: a lactone-fused benzene ring
        ("Coumarin", "O=C1Oc2ccccc2C1"),
        # Flavone motif (common in flavonoids)
        ("Flavone", "c1cc2oc(=O)cc(c2c1)"),
        # Isoflavonoid motif:
        ("Isoflavonoid", "c1ccc2c(c1)cc(=O)oc2"),
        # Generic unsaturated C6-C3: benzene ring attached to a C=C unit then any carbon.
        ("Unsaturated C6-C3", "c1ccccc1/C=C/[CH2]")
    ]

    # Try each SMARTS pattern.
    for name, smarts in patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # Skip if the SMARTS pattern is not valid.
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches the {name} pattern, likely a phenylpropanoid."

    # Fallback: if no specialized pattern matched, require the presence
    # of a benzene ring. (This is a weak test but helps avoid false positives
    # for non-aromatic compounds.)
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene):
        return False, "No benzene (aromatic) ring found, unlikely to be a phenylpropanoid."

    return False, "No defined phenylpropanoid substructure pattern matched."

# Example usage (for testing; can be commented out in production):
if __name__ == "__main__":
    test_smiles = [
        "CCOC(=O)\\C=C\\c1ccccc1",                 # Ethyl cinnamate: Cinnamoyl Ester
        "COC(=O)\\C=C\\c1ccc(OC)cc1",               # methyl 4-methoxycinnamate: Cinnamoyl Ester
        "Oc1cc(=O)oc2ccccc12",                      # 4-hydroxycoumarin: Coumarin pattern
        "C1=C(C(=C2C(=C1)C=CC(O2)=O)CC=C(C)C)O",     # osthenol: may be caught by unsaturated C6-C3
        "O1[C@@H]([C@H]([C@@H]([C@@H]1C2=CC(OC)=C(O)C=C2)C)C)C3=CC=4OCOC4C=C3",  # Austrobailignan 7: flavonoid-type
        "COc1cc(C=CC)ccc1O"                         # isoeugenol: matches unsaturated C6-C3
    ]

    for smi in test_smiles:
        result, reason = is_phenylpropanoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")