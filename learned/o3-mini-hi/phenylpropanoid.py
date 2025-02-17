"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: Phenylpropanoid – any aromatic compound based on a phenylpropane (C6–C3) skeleton.

This heuristic classifier uses several substructure SMARTS patterns to detect:
  • A cinnamoyl motif (both forward and reverse orientations),
  • The coumarin lactone ring,
  • A flavone (and isoflavonoid) fused-ring core,
  • And a generic C6–C3 chain (both saturated and unsaturated).

Phenylpropanoids are very diverse so this classifier is not perfect.
If none of the defined patterns match we finally verify that a benzene ring is at least present.
"""

from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    
    Uses a set of SMARTS patterns that capture substructures such as a cinnamoyl unit
    (which may be in either orientation), coumarin, flavone, isoflavonoid,
    or a benzene directly attached to a three-carbon chain.
    
    Args:
        smiles (str): A SMILES string representing the molecule.
    
    Returns:
        bool: True if classified as a phenylpropanoid, False otherwise.
        str: Explanation for the classification decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a list of (pattern name, SMARTS) tuples.
    # Note: The cinnamoyl patterns are defined in both directions.
    patterns = [
        # Cinnamoyl patterns: benzene attached via a trans double bond to a carbonyl as acid or ester.
        ("Cinnamoyl Ester", "c1ccccc1/C=C/C(=O)O"),
        ("Cinnamoyl Acid", "O=C(O)/C=C/c1ccccc1"),
        # Coumarin (lactone fused benzene ring) pattern:
        ("Coumarin", "O=C1Oc2ccccc2C1"),
        # Flavone motif (common in flavonoids)
        ("Flavone", "c1cc2oc(=O)cc(c2c1)"),
        # Isoflavonoid motif:
        ("Isoflavonoid", "c1ccc2c(c1)cc(=O)oc2"),
        # Generic unsaturated C6-C3: benzene attached to a C=C and then any carbon.
        ("Unsaturated C6-C3", "c1ccccc1/C=C/[C]"),
        # Generic saturated C6-C3: benzene attached to exactly three aliphatic carbons.
        ("Saturated C6-C3", "c1ccccc1CCC")
    ]
    
    # Try each SMARTS pattern.
    for name, smarts in patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # skip if the SMARTS pattern is invalid
        if mol.HasSubstructMatch(pattern):
            # If a generic pattern, we can be less strict.
            return True, f"Matches the {name} pattern, likely a phenylpropanoid."
            
    # Fallback: Require that at least one benzene ring is present.
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene):
        return False, "No benzene (aromatic) ring found, unlikely to be a phenylpropanoid."
    
    return False, "No common phenylpropanoid substructure pattern matched."

# Example usage (for testing – can be removed or commented out in production):
if __name__ == "__main__":
    test_smiles = [
        "CCOC(=O)\\C=C\\c1ccccc1",               # Ethyl cinnamate: should match Cinnamoyl patterns.
        "COC(=O)\\C=C\\c1ccc(OC)cc1",             # methyl 4-methoxycinnamate: cinnamoyl.
        "Oc1cc(=O)oc2ccccc12",                    # 4-hydroxycoumarin: coumarin pattern.
        "C1=C(C(=C2C(=C1)C=CC(O2)=O)CC=C(C)C)O",   # osthenol (fused ring): may match flavone-type.
        "O1C([C@H]([C@@H]([C@@H]([C@H]1C2=CC(OC)=C(O)C=C2)C)C)C)C3=CC=4OCOC4C=C3"  # Austrobailignan 7.
    ]
    for smi in test_smiles:
        result, reason = is_phenylpropanoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")