"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: Flavonoid
Definition: A flavonoid is any member of the 'superclass' flavonoids whose 
skeleton is based on 1-benzopyran with an aryl substituent at position 2.
This heuristic implementation checks for three classical substructures:
  1. Flavone: 2-phenylchromen-4-one scaffold.
  2. Flavanone: A related scaffold with a saturated central ring.
  3. Flavanol: A variant with an –OH substituent on the heterocycle.
Note: Many flavonoids bear additional substituents (e.g. sugars) but if 
the core is present, the molecule is classified as a flavonoid.
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is (putatively) a flavonoid based on its SMILES string.
    
    This function uses substructure searches with SMARTS queries that capture
    classical flavonoid cores (e.g. flavone, flavanone, flavanol). These patterns
    have an internal benzopyran system with an exocyclic aromatic ring (at position 2).
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule appears to contain a flavonoid scaffold, else False.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for typical flavonoid cores.
    # 1. Flavone: 2-phenylchromen-4-one. Example: quercetin aglycone.
    pattern_flavone = Chem.MolFromSmarts("c1ccc(cc1)-c2oc(=O)c3ccccc3c2")
    
    # 2. Flavanone: Similar to flavone but the pyran ring is saturated except for the carbonyl.
    # This pattern represents a 2-phenylchroman-4-one.
    pattern_flavanone = Chem.MolFromSmarts("c1ccc(cc1)-C2CC(=O)c3ccccc3O2")
    
    # 3. Flavanol: Like flavanone but with a hydroxyl instead of carbonyl at position 4.
    # Example: catechin.
    pattern_flavanol = Chem.MolFromSmarts("c1ccc(cc1)-C2CC(O)c3ccccc3O2")
    
    # Try matching each pattern; if one is found, we classify as flavonoid.
    if mol.HasSubstructMatch(pattern_flavone):
        return True, "Molecule contains a flavone (2-phenylchromen-4-one) flavonoid scaffold."
    if mol.HasSubstructMatch(pattern_flavanone):
        return True, "Molecule contains a flavanone (2-phenylchroman-4-one) flavonoid scaffold."
    if mol.HasSubstructMatch(pattern_flavanol):
        return True, "Molecule contains a flavanol (2-phenylchroman-4-ol) flavonoid scaffold."
    
    # If none of the patterns matches, we return False.
    return False, "No recognizable flavonoid scaffold (based on classical flavone/flavanone/flavanol cores) was detected."

# Example usage (for testing)
if __name__ == "__main__":
    # (R)-naringenin is a prototypical flavanone.
    test_smiles = "Oc1ccc(cc1)[C@H]1CC(=O)c2c(O)cc(O)cc2O1"
    result, reason = is_flavonoid(test_smiles)
    print(result, ":", reason)