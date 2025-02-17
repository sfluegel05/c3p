"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: Limonoid class (tetranortriterpenoids that are highly oxygenated and contain or derive from
a 4,4,8-trimethyl-17-furanylsteroid skeleton). Note: This is a heuristic filter.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    
    Heuristics used:
      - The molecule must be parseable.
      - It should be polycyclic (at least 4 rings) to resemble the steroid/terpenoid skeleton.
      - The carbon count should fall in a range suggestive of a triterpenoid (roughly 25-35 carbons)
        since tetranortriterpenoids have lost 4 carbons from a 30-carbon precursor.
      - The molecule should be highly oxygenated; here we require that the oxygen/carbon ratio > 0.2.
      - It should contain a furan ring (SMARTS: 5-membered aromatic ring containing one O).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as limonoid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Count the number of rings using SSSR (Smallest Set of Smallest Rings)
    num_rings = Chem.GetSSSR(mol)
    if num_rings < 4:
        return False, f"Insufficient rings found (only {num_rings} rings; expected at least 4 for a terpenoid skeleton)."
    
    # Count carbon and oxygen atoms (atomic numbers 6 and 8, respectively)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check if carbon count is in the plausible range for triterpenoids/tetranortriterpenoids.
    if not (25 <= carbon_count <= 35):
        return False, f"Carbon count ({carbon_count}) not within the expected range for a triterpenoid/tetranortriterpenoid."
    
    # Check for high oxygenation: using an arbitrary oxygen:carbon ratio threshold.
    # Limonoids are known to be highly oxygenated.
    ratio = oxygen_count / carbon_count if carbon_count > 0 else 0
    if ratio < 0.2:
        return False, f"Oxygenation level low (O/C ratio = {ratio:.2f}; expected > 0.2 for highly oxygenated molecules)."
    
    # Look for a furan ring substructure.
    # Furan: aromatic 5-membered ring containing one oxygen (pattern: "c1occc1").
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found, which is expected in a 4,4,8-trimethyl-17-furanylsteroid skeleton or derivative."
    
    # If all conditions are met then the molecule qualifies (heuristically) as a limonoid.
    return True, "Meets heuristic criteria for a limonoid (polycyclic, appropriate carbon count, highly oxygenated, and contains a furan ring)."

# Optional: Testing the function on an example SMILES (one of the provided limonoid examples)
if __name__ == "__main__":
    # Example: Deacetylnomilin (one of the provided limonoid examples)
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_limonoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)