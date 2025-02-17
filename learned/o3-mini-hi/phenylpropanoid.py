"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: Phenylpropanoid – any aromatic compound based on a phenylpropane (C6-C3) skeleton.

In this revised approach we use several SMARTS patterns:
  • Cinnamoyl Ester and Acid motifs: common in cinnamate esters/acids.
  • Coumarin: a benzopyran-2-one motif (with a relaxed pattern to allow substituents).
  • Flavone and Isoflavonoid: for flavonoid cores.
  • Allylbenzene: a benzene ring directly attached to a C=CC fragment
which helps catch phenylpropanoids lacking a carbonyl group (e.g. isoeugenol).
  
If none of these patterns match but the molecule does not even have a benzene ring,
we label it as non-phenylpropanoid.
"""

from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.

    Uses a set of SMARTS patterns that capture key substructures typical of phenylpropanoid subclasses:
      • Cinnamoyl Acid/Ester motifs,
      • Coumarin (benzopyran-2-one),
      • Flavone and Isoflavonoid cores,
      • And a benzene ring with an allylic (C=CC) side chain.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phenylpropanoid, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a list of SMARTS patterns (name, SMARTS string)
    patterns = [
        # Cinnamoyl ester: benzene ring attached via a trans double bond to a carbonyl (ester)
        ("Cinnamoyl Ester", "c1ccc(/C=C/C(=O)O)cc1"),
        # Cinnamoyl acid: similar but as an acid
        ("Cinnamoyl Acid", "O=C(O)/C=C/c1ccccc1"),
        # Coumarin: benzopyran-2-one; allow non‐critical substituents by using a relaxed pattern
        ("Coumarin", "O=c1ccc2occc2c1"),
        # Flavone motif (common in flavonoids)
        ("Flavone", "c1cc2oc(=O)cc(c2c1)"),
        # Isoflavonoid motif:
        ("Isoflavonoid", "c1ccc2c(c1)cc(=O)oc2"),
        # Allylbenzene substructure: benzene ring with a C=CC side chain
        ("Allylbenzene", "c1ccc(C=CC)cc1")
    ]

    # Test for at least one high-confidence phenylpropanoid pattern.
    for name, smarts in patterns:
        substruct = Chem.MolFromSmarts(smarts)
        if substruct is None:
            continue  # skip invalid SMARTS
        if mol.HasSubstructMatch(substruct):
            return True, f"Matches the {name} pattern, likely a phenylpropanoid."

    # Fallback: if no specialized pattern matched, require the presence of a benzene ring.
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene):
        return False, "No benzene (aromatic) ring found, unlikely to be a phenylpropanoid."

    return False, "No defined phenylpropanoid substructure pattern matched."

# Example usage for testing:
if __name__ == "__main__":
    test_smiles = [
        "CCOC(=O)\\C=C\\c1ccccc1",                 # Ethyl cinnamate (Cinnamoyl Ester)
        "COC(=O)\\C=C\\c1ccc(OC)cc1",               # methyl 4-methoxycinnamate (Cinnamoyl Ester)
        "Oc1cc(=O)oc2ccccc12",                      # 4-hydroxycoumarin (Coumarin)
        "C1=C(C(=C2C(=C1)C=CC(O2)=O)CC=C(C)C)O",     # osthenol (should be caught by a flavonoid/coumarin-type context)
        "O1[C@@H]([C@H]([C@@H]([C@@H]1C2=CC(OC)=C(O)C=C2)C)C)C3=CC=4OCOC4C=C3",  # Austrobailignan 7 (flavonoid core)
        "COc1cc(C=CC)ccc1O"                         # isoeugenol (Allylbenzene)
    ]

    for smi in test_smiles:
        result, reason = is_phenylpropanoid(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("")