"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: Phenylpropanoid – any aromatic compound based on a phenylpropane (C6–C3) skeleton.
Improved heuristics:
  • Look for a benzene ring attached to exactly a three‐carbon chain (the C6–C3 skeleton)
    either as a saturated chain or one containing a C=C bond (cinnamoyl-type).
  • Look for the lactone/fused ring pattern of coumarins.
  • Look for the core patterns of flavonoids and isoflavonoids.
A secondary check ensures that a benzene ring is present.
Note: Due to the diversity of phenylpropanoid derivatives, this is a heuristic classifier.
"""

from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    The heuristics try to match:
      - A benzene ring directly attached to a three‐carbon chain (saturated or with a C=C double bond);
      - A cinnamoyl-type substructure (where the propanoid unit terminates with a carbonyl);
      - The coumarin lactone fused ring system;
      - A flavonoid or isoflavonoid fused ring core.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a phenylpropanoid, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # List of tuples: (pattern name, SMARTS)
    smarts_patterns = [
        # Saturated C6-C3: benzene ring attached to exactly three methylene groups.
        ("Saturated C6-C3", "c1ccccc1-[CH2]-[CH2]-[CH2]"),
        # Unsaturated C6-C3: benzene ring attached via a single bond to a C=C and then any carbon.
        ("Unsaturated C6-C3", "c1ccccc1/C=C/[C]"),
        # Cinnamoyl-type: benzene ring attached via a C=C to a carbonyl carbon.
        ("Cinnamoyl", "c1ccccc1/C=C/C(=O)"),
        # Coumarin motif:
        ("Coumarin", "O=C1Oc2ccccc2C1"),
        # Flavonoid core:
        ("Flavonoid", "c1ccc2c(c1)oc(=O)c3ccccc23"),
        # Isoflavonoid core:
        ("Isoflavonoid", "c1ccc2c(c1)cc(=O)oc2")
    ]
    
    # Iterate through our defined substructure patterns.
    for name, smarts in smarts_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # if SMARTS parsing fails, skip pattern
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            # For the generic C6-C3 chain patterns, perform an additional check to see if the attachment is “clean”
            if name in ["Saturated C6-C3", "Unsaturated C6-C3", "Cinnamoyl"]:
                # We expect that the bond connecting the benzene and the chain originates
                # from a benzene carbon that is not overly substituted.
                # For each match, the first atom in the substructure (index from our SMARTS) should be
                # the benzene carbon. Its degree (ignoring hydrogens) should be 3 (two neighbours in the ring plus one chain bond).
                valid_match = False
                for match in matches:
                    # In our SMARTS, atom 0 is the benzene atom (part of c1ccccc1).
                    benzene_atom = mol.GetAtomWithIdx(match[0])
                    heavy_deg = sum(1 for nbr in benzene_atom.GetNeighbors() if nbr.GetAtomicNum() > 1)
                    if heavy_deg == 3:
                        valid_match = True
                        break
                if not valid_match:
                    continue  # Skip this pattern if none of the matches pass the attachment check.
            return True, f"Matches the {name} pattern, likely a phenylpropanoid."
    
    # As a fallback, check for at least one benzene ring.
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene):
        return False, "No benzene (aromatic) ring found, unlikely to be a phenylpropanoid."
    
    return False, "No common phenylpropanoid substructure pattern matched."

# Example usage (can be removed in production):
if __name__ == "__main__":
    test_smiles = [
        "CCOC(=O)\\C=C\\c1ccccc1",   # Ethyl cinnamate – should match Unsaturated C6-C3 or Cinnamoyl
        "COC(=O)\\C=C\\c1ccc(OC)cc1", # Methyl 4-methoxycinnamate – unsaturated C6-C3/cinnamoyl
        "Oc1cc(=O)oc2ccccc12",        # 4-hydroxycoumarin – coumarin motif
        "c1=C(c(c2C(=C1)C=CC(O2)=O)CC=C(C)C)O",  # osthenol – expected to match one of the fused ring patterns
    ]
    
    for smi in test_smiles:
        result, reason = is_phenylpropanoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")