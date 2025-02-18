"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: Any 1-benzopyran with an aryl substituent at position 3 (isoflavonoid).

An isoflavonoid is defined as a benzopyran system (a chromene)
where a fused 6-membered aromatic ring (benzene) is joined to a 6-membered ring containing one oxygen
and an exocyclic aryl substituent is attached at the position (position 3) of the pyran ring.
This implementation uses two SMARTS queries – one for an sp2 variant (isoflavone-like)
and one for an sp3 variant (isoflavan-like) – to better capture the range of isoflavonoid scaffolds.
"""

from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.
    This implementation uses two SMARTS substructure patterns:
      • Pattern 1: "c1ccc2c(c1)oc(c2)c3ccccc3" 
        => represents a chromene core (benzopyran) in which the carbon at position 3 remains aromatic 
           (as in many isoflavones).
      • Pattern 2: "c1ccc2c(c1)oc(c2)C(c3ccccc3)" 
        => represents the same core but with an aliphatic (sp3) carbon at the 3‑position (as in isoflavans).
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an isoflavonoid, False otherwise.
        str: Explanation of the classification.
    """
    
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compile two SMARTS patterns to detect the isoflavonoid core.
    # Pattern 1: The 3-position remains aromatic (isoflavone-type).
    smarts1 = "c1ccc2c(c1)oc(c2)c3ccccc3"
    pattern1 = Chem.MolFromSmarts(smarts1)
    
    # Pattern 2: The 3-position is aliphatic (typical of isoflavan-type scaffolds).
    smarts2 = "c1ccc2c(c1)oc(c2)C(c3ccccc3)"
    pattern2 = Chem.MolFromSmarts(smarts2)
    
    # Check if the molecule matches either of the patterns.
    if mol.HasSubstructMatch(pattern1) or mol.HasSubstructMatch(pattern2):
        return True, ("Molecule contains a fused benzopyran (chromene) core with an exocyclic aryl "
                      "substituent at the expected 3‐position, consistent with an isoflavonoid scaffold.")
    
    # If neither SMARTS pattern matches, then we do not classify it as an isoflavonoid.
    return False, "Scaffold not recognized as isoflavonoid (missing benzopyran core with 3-aryl substituent)."

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        # Some examples — note that these are simplified or may be fragments.
        "c1ccc2c(c1)oc(c2)c3ccccc3",  # Expected: isoflavonoid (pattern 1)
        "c1ccc2c(c1)oc(c2)C(c3ccccc3)",  # Expected: isoflavonoid (pattern 2)
        "O=C1C(=O)C2=CC=CC=C2O1",  # Not isoflavonoid
    ]
    for smi in test_smiles:
        res, reason = is_isoflavonoid(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")