"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
#!/usr/bin/env python
"""
Classifies: Phosphatidylcholine (A glycerophosphocholine that is glycero-3-phosphocholine bearing two acyl substituents at positions 1 and 2.)
"""

from rdkit import Chem

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    The logic first attempts to match a SMARTS pattern intended to capture a glycerol backbone
    bearing two acyl chains (via ester bonds) and a phosphocholine headgroup.
    If that fails, it falls back to checking for a choline fragment, a phosphate group,
    and at least two acyl ester groups.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as phosphatidylcholine, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Primary approach: Use a SMARTS pattern that captures the expected phosphatidylcholine backbone.
    # This pattern looks for a central carbon with two esterified substituents and a phosphate
    # ester connected to a choline group.
    # The pattern string is: 
    #   C(OC(=O)[*])(OC(=O)[*])OP(=O)([O-])OCC[N+](C)(C)C
    # where [*] is a wildcard for the acyl chains.
    pc_pattern = "C(OC(=O)[*])(OC(=O)[*])OP(=O)([O-])OCC[N+](C)(C)C"
    pc_mol = Chem.MolFromSmarts(pc_pattern)
    if pc_mol is not None and mol.HasSubstructMatch(pc_mol):
        return True, "Matches phosphatidylcholine backbone SMARTS pattern."

    # Fallback approach:
    # (1) Check for the choline fragment (OCC[N+](C)(C)C)
    choline_pattern = Chem.MolFromSmarts("OCC[N+](C)(C)C")
    if choline_pattern is None or not mol.HasSubstructMatch(choline_pattern):
        return False, "Missing choline head group (OCC[N+](C)(C)C)."

    # (2) Check for a phosphate group - require at least one phosphorus atom.
    phosphate_pattern = Chem.MolFromSmarts("[P]")
    if phosphate_pattern is None or not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing phosphate group."

    # (3) Check for at least two acyl ester bonds through the fragment "OC(=O)"
    acyl_pattern = Chem.MolFromSmarts("OC(=O)")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 2:
        return False, f"Found {len(acyl_matches)} acyl ester group(s); at least two are required for phosphatidylcholine."

    # If all fallback conditions are met, assume the molecule is phosphatidylcholine.
    return True, "Contains choline head group, phosphate and at least two acyl ester groups consistent with phosphatidylcholine structure."

# Example usage:
if __name__ == '__main__':
    # Test using one of the provided PC SMILES strings.
    test_smiles = "P(OC[C@@H](COC(CCCCCCCCCCCCCCCCC)=O)OC(=O)CCCC)(=O)(OCC[N+](C)(C)C)[O-]"
    result, explanation = is_phosphatidylcholine(test_smiles)
    print("Result:", result)
    print("Explanation:", explanation)