"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: Flavones
A flavone is any flavonoid having a 2-aryl-1-benzopyran-4-one (2-phenylchromen-4-one) skeleton.
This function determines if a given SMILES string contains that scaffold – even if substitutions (such as sugars)
are present.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone (or a substituted flavone derivative) based on its SMILES string.
    
    The criteria are simplified:
      1. The molecule must contain a 2-phenylchromen-4-one (flavone) core as detected by a SMARTS pattern.
      2. We allow extra substituents (for example sugars) and do not require exactly 15 heavy atoms in the match.
      3. As an extra check, we compute the Murcko scaffold and ensure that the core ring system 
         consists of at least 2–3 rings, which is typical for flavonoids.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as having a flavone core, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the 2-aryl-1-benzopyran-4-one (flavone) core.
    # This pattern is intentionally written in a relaxed way: it will match the core even if extra substituents are present.
    flavone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)c(c2)-c3ccccc3")
    if flavone_pattern is None:
        return False, "Error in SMARTS definition for flavone core"
    
    # Check if the molecule contains the flavone core as a substructure.
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "Does not contain a 2-phenylchromen-4-one (flavone) core"
    
    # Optionally, examine the Murcko scaffold to see if the backbone is consistent with flavonoids.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    ring_info = scaffold.GetRingInfo()
    n_rings = ring_info.NumRings() if ring_info is not None else 0
    # Typically, a flavone has (at least) two fused rings and a third aromatic ring (the 2-phenyl substituent).
    if n_rings < 2:
        return False, f"Scaffold only has {n_rings} ring(s), which seems insufficient for a flavone core"
    
    # If we made it here, we assume that the flavone core is present.
    return True, "Contains a flavone (2-phenylchromen-4-one) core (allowing for substituted derivatives)"

# For demonstration purposes only:
if __name__ == "__main__":
    # Test with a known flavone (3-methoxyapigenin) and print results.
    test_smiles = "COc1c(oc2cc(O)cc(O)c2c1=O)-c1ccc(O)cc1"
    result, reason = is_flavones(test_smiles)
    print("Test molecule:", test_smiles)
    print("Is flavone?", result)
    print("Reason:", reason)