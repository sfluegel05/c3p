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
      - Must be a parseable molecule.
      - Should be polycyclic (at least 4 rings) to resemble a typical steroid/terpenoid skeleton.
      - The carbon count should be within a plausible range for tetranortriterpenoids (roughly 25-35 carbons).
      - The molecule should be highly oxygenated; here we require an oxygen/carbon ratio > 0.2.
      - It should contain a furan ring (SMARTS: 5-membered aromatic ring containing one oxygen), as part
        of the 4,4,8-trimethyl-17-furanylsteroid skeleton or a derivative thereof.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the heuristic criteria for a limonoid, False otherwise.
        str: A message describing the reasoning for classification.
    """
    # Attempt to create an RDKit Mol object from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Calculate the number of rings using a reliable RDKit descriptor.
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 4:
        return False, f"Insufficient rings found (only {num_rings} rings; expected at least 4 for a terpenoid skeleton)."
    
    # Count the number of carbon and oxygen atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check if the carbon count is within a plausible range for tetranortriterpenoids.
    if not (25 <= carbon_count <= 35):
        return False, f"Carbon count ({carbon_count}) not within the expected range (25-35) for a tetranortriterpenoid."
    
    # Calculate the oxygen to carbon ratio.
    ratio = oxygen_count / carbon_count if carbon_count > 0 else 0
    if ratio < 0.2:
        return False, f"Oxygenation level low (O/C ratio = {ratio:.2f}; expected > 0.2 for highly oxygenated molecules)."
    
    # Look for a furan ring substructure.
    # Furan: aromatic 5-membered ring containing one oxygen; SMARTS pattern "c1occc1".
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found, which is expected in a 4,4,8-trimethyl-17-furanylsteroid skeleton or derivative."
    
    return True, ("Meets heuristic criteria for a limonoid: polycyclic structure, appropriate carbon count, high oxygenation, "
                  "and contains a furan ring.")

# Optional testing when running as a script:
if __name__ == "__main__":
    # Test with one of the provided example SMILES (Deacetylnomilin)
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_limonoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)