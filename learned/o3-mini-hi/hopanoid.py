"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: A hopanoid is defined as a triterpenoid based on a hopane skeleton.
We use the molecule’s Murcko scaffold as a proxy for the fused polycyclic core.
The new approach checks that:
  - The scaffold has a carbon count between 20 and 34 (hopane is C30,
    but derivative scaffolds may lose some carbons).
  - The scaffold contains at least 5 rings.
  - The scaffold is “carbon pure” (i.e. very few heteroatoms).
  - The overall molecule has enough carbons to be a triterpenoid.
This heuristic is not perfect but attempts to reduce both false negatives and false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_hopanoid(smiles: str):
    """
    Determines whether a molecule is a hopanoid based on its SMILES string.
    A hopanoid is defined as a triterpenoid based on a hopane skeleton.
    The algorithm:
      1. Extracts the Murcko scaffold.
      2. Calculates the number of carbon atoms in the scaffold.
      3. Calculates the ratio of non-carbon atoms in the scaffold.
      4. Counts the number of rings in the scaffold.
      5. Performs an overall carbon-count sanity check.
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule is classified as a hopanoid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Get the Murcko scaffold (the core ring system).
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Could not extract a Murcko scaffold."

    # Count the number of carbon atoms in the scaffold.
    scaffold_c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetSymbol() == "C")
    # Count the total number of atoms in the scaffold.
    scaffold_total = scaffold.GetNumAtoms()
    # Calculate fraction of atoms that are NOT carbon.
    non_c_count = scaffold_total - scaffold_c_count
    hetero_ratio = non_c_count / scaffold_total if scaffold_total > 0 else 1.0

    # Get the number of rings in the scaffold.
    ring_info = scaffold.GetRingInfo()
    num_rings = ring_info.NumRings() if ring_info is not None else 0

    # Count total number of carbons in the whole molecule.
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")

    # Heuristic checks:
    # (1) Accept scaffold carbon count between 20 and 34.
    if not (20 <= scaffold_c_count <= 34):
        return False, (f"Murcko scaffold has {scaffold_c_count} carbons; "
                       "expected between 20 and 34 for a hopane skeleton.")
    
    # (2) The scaffold should have at least 5 rings, reflecting the pentacyclic core.
    if num_rings < 5:
        return False, f"Murcko scaffold has {num_rings} rings; at least 5 rings required."
    
    # (3) The scaffold should be mostly carbon (i.e. low heteroatom ratio).
    # Hopane cores are almost entirely carbon. Allow a small margin (e.g. 15%).
    if hetero_ratio > 0.15:
        return False, (f"Murcko scaffold contains {non_c_count} non-carbon atoms out of {scaffold_total} "
                       f"({hetero_ratio*100:.1f}%); expected a mostly carbon scaffold.")
    
    # (4) Overall, the molecule should have enough carbons to be a triterpenoid.
    if total_c < 25:
        return False, "Total carbon count is too low for a triterpenoid."
    
    return True, ("Molecule has a fused polycyclic scaffold with a mostly carbon backbone "
                  f"({scaffold_c_count} C, {num_rings} rings), consistent with a hopane-based triterpenoid.")

# Example usage:
if __name__ == "__main__":
    # Test with an example hopanoid (e.g., dustanin):
    dustanin = "CC(C)(O)[C@H]1CC[C@@]2(C)[C@H]1C[C@H](O)[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3CC[C@@]12C"
    result, reason = is_hopanoid(dustanin)
    print("Is hopanoid:", result)
    print("Reason:", reason)