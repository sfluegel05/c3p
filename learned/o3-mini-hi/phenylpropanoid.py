"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: Phenylpropanoid – any aromatic compound based on a phenylpropane (C6–C3) skeleton.
This heuristic classification checks if the input molecule contains typical substructures
such as a benzene ring linked to a three-carbon chain (saturated or with unsaturation) or common
fused ring systems found in coumarins, flavonoids, and isoflavonoids.
"""

from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    The heuristic is based on the presence of a benzene ring with an attached three-carbon moiety
    (the C6–C3 skeleton) or characteristic fused ring systems in coumarins/flavonoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a phenylpropanoid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # List of tuples (pattern_name, SMARTS) for typical phenylpropanoid substructures.
    # Note: These are heuristic patterns that capture common C6-C3 and fused ring motifs.
    smarts_patterns = [
        # Saturated phenylpropanoid: benzene ring attached to a three-carbon chain (saturated)
        ("Saturated C6-C3", "c1ccccc1CCC"),
        # Unsaturated (cinnamoyl) type: benzene ring attached to a C=C chain then a carbon (e.g., cinnamates)
        ("Unsaturated C6-C3", "c1ccccc1C=CC"),
        # Coumarin core (e.g., 4-hydroxycoumarin); note that many coumarins have this characteristic pattern.
        ("Coumarin", "Oc1cc(=O)oc2ccccc12"),
        # Flavonoid core pattern (a chromanone fused with benzene) – typical for flavonoids.
        ("Flavonoid", "c1ccc2c(c1)oc(=O)c3ccccc23"),
        # Isoflavonoid – a variant present in many natural phenylpropanoids.
        ("Isoflavonoid", "c1ccc2c(c1)cc(=O)oc2")
    ]
    
    # Iterate over the SMARTS patterns to check if any is found in the molecule
    for name, smarts in smarts_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        # If the molecule has the substructure, we classify it as a phenylpropanoid.
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return True, f"Matches the {name} pattern, likely a phenylpropanoid."
    
    # As an extra sanity check, require the molecule contains at least one benzene ring.
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene):
        return False, "No benzene (aromatic) ring found, unlikely to be a phenylpropanoid."
    
    # If none of the typical C6-C3 or fused ring systems match, then return false.
    return False, "No common phenylpropanoid substructure pattern matched."

# Example usage:
if __name__ == "__main__":
    smiles_examples = [
        "CCOC(=O)\\C=C\\c1ccccc1",  # Ethyl cinnamate (should match Unsaturated C6-C3)
        "COC(=O)\\C=C\\c1ccc(OC)cc1",  # Methyl 4-methoxycinnamate
        "Oc1cc(=O)oc2ccccc12"  # 4-hydroxycoumarin
    ]
    
    for smi in smiles_examples:
        result, reason = is_phenylpropanoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")