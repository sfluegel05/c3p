"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: Flavones
A flavone is any flavonoid having a 2-aryl-1-benzopyran-4-one (2-phenylchromen-4-one) skeleton.
This version improves on the earlier approach by stripping away peripheral substituents using the Murcko scaffold,
and then checking the scaffold with a stricter SMARTS pattern.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone derivative based on its SMILES string.
    
    The new strategy is:
      1. Parse the molecule.
      2. Compute its Murcko scaffold to remove common substituents (such as sugars) that might mask the core.
      3. Use a stricter SMARTS pattern that matches a 2-phenylchromen-4-one.
         The pattern "c1ccccc1-C2=CC(=O)Oc3ccccc23" ensures a fused benzopyran-4-one core 
         (rings A and C) with a phenyl substituent (ring B) attached at C2.
      4. As an extra safeguard, we verify that the extracted scaffold has at least 3 rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a flavone derivative, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute the Murcko scaffold to extract the core structure.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        # Fallback: use the original molecule if scaffold extraction fails.
        scaffold = mol
        
    # Count rings in the scaffold.
    ring_info = scaffold.GetRingInfo()
    n_rings = ring_info.NumRings() if ring_info is not None else 0
    if n_rings < 3:
        return False, f"Scaffold only has {n_rings} ring(s), which is insufficient for a flavone core"
    
    # Define a stricter SMARTS pattern for a 2-phenylchromen-4-one scaffold.
    # This explicitly requires a phenyl ring attached to a benzopyran-4-one system.
    flavone_pattern = Chem.MolFromSmarts("c1ccccc1-C2=CC(=O)Oc3ccccc23")
    if flavone_pattern is None:
        return False, "Error in SMARTS definition for flavone core"
    
    # Check if the scaffold contains the flavone core.
    if not scaffold.HasSubstructMatch(flavone_pattern):
        return False, "Does not contain a 2-phenylchromen-4-one (flavone) core in its scaffold"
    
    # Passed all tests: the molecule appears to be a flavone (or substituted flavone derivative).
    return True, "Scaffold contains a 2-phenylchromen-4-one core typical of flavones (after removal of substituents)"

# For demonstration purposes only:
if __name__ == "__main__":
    # A few test examples:
    test_examples = [
        ("COc1c(oc2cc(O)cc(O)c2c1=O)-c1ccc(O)cc1", "3-methoxyapigenin"),
        ("OC[C@@H]1O[C@H](Oc2c(oc3cc(O)cc(O)c3c2=O)-c2cc(O)c(O)c(O)c2)[C@@H](O)[C@H](O)[C@@H]1O", "myricetin 3-O-beta-L-galactopyranoside")
    ]
    
    for smiles, name in test_examples:
        result, reason = is_flavones(smiles)
        print(f"Test molecule: {name}")
        print(f"SMILES: {smiles}")
        print(f"Is flavone? {result}")
        print(f"Reason: {reason}\n")