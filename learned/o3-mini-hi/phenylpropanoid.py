"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: Phenylpropanoid – any aromatic compound based on a phenylpropane (C6–C3) skeleton.
Heuristics:
  • Look for a benzene ring directly attached to a three‐carbon side chain using only CH2 atoms
    (the saturated C6–C3 skeleton).
  • Look for an unsaturated form where at least one C–C bond in the chain is a double bond (cinnamoyl-type).
  • Look for the lactone/fused ring pattern of coumarins.
  • Look for the common motifs of flavonoids and isoflavonoids.
Keep in mind that many phenylpropanoid derivatives may require additional context.
"""

from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    The heuristic checks for a benzene ring with an attached three-carbon chain (the C6–C3 skeleton)
    and/or characteristic fused ring systems found in coumarins, flavonoids and isoflavonoids.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a phenylpropanoid, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # List of tuples (pattern_name, SMARTS) for typical phenylpropanoid substructures.
    # We now attempt to limit the match by requiring that the three-carbon chain consists solely of CH2
    # atoms (for the saturated case) and one C=C bond for the unsaturated case.
    smarts_patterns = [
        # Saturated phenylpropanoid: benzene ring attached to a straight chain of three CH2 groups.
        ("Saturated C6-C3", "c1ccccc1[CH2][CH2][CH2]"),
        # Unsaturated (cinnamoyl-type): benzene ring attached to a chain with one double bond.
        ("Unsaturated C6-C3", "c1ccccc1[CH2]C=CC"),
        # Coumarin-like pattern: lactone-fused benzene ring.
        ("Coumarin", "O=C1Oc2ccccc2C1"),
        # Flavonoid core: a benzopyranone fused with a benzene ring.
        ("Flavonoid", "c1ccc2c(c1)oc(=O)c3ccccc23"),
        # Isoflavonoid: a relocation of the attachment; common in many natural products.
        ("Isoflavonoid", "c1ccc2c(c1)cc(=O)oc2")
    ]
    
    # Iterate over the defined SMARTS patterns
    for name, smarts in smarts_patterns:
        substruct = Chem.MolFromSmarts(smarts)
        if substruct is None:
            continue  # Skip if SMARTS is not parsed correctly.
        if mol.HasSubstructMatch(substruct):
            return True, f"Matches the {name} pattern, likely a phenylpropanoid."
    
    # As a secondary check, we require at least one benzene ring.
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene):
        return False, "No benzene (aromatic) ring found, unlikely to be a phenylpropanoid."
    
    # No typical phenylpropanoid substructure was found.
    return False, "No common phenylpropanoid substructure pattern matched."

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "CCOC(=O)\\C=C\\c1ccccc1", # Ethyl cinnamate – unsaturated C6-C3
        "COC(=O)\\C=C\\c1ccc(OC)cc1", # Methyl 4-methoxycinnamate – unsaturated C6-C3
        "Oc1cc(=O)oc2ccccc12",       # 4-hydroxycoumarin – coumarin motif
        "c1=C(c(c2C(=C1)C=CC(O2)=O)CC=C(C)C)O", # osthenol – might match one of the fused ring patterns
    ]
    
    for smi in test_smiles:
        result, reason = is_phenylpropanoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")