"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: A hopanoid is defined as a triterpenoid based on a hopane skeleton.
This program uses a Murcko scaffold to try to capture a fused polycyclic core 
that should have roughly 30 carbons (hopane has C30) and at least 5 rings.
Note: This is a heuristic approach and may not classify every hopanoid correctly.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is defined as a triterpenoid based on a hopane skeleton.
    The approach here extracts the molecule's Murcko scaffold, then checks:
      - The number of carbons in the scaffold (should be roughly that of the hopane core, ~30)
      - That the scaffold contains at least 5 rings (consistent with a pentacyclic hopane skeleton)
      - A sanity check that the overall molecule contains enough carbons 
        (suggesting it is indeed a triterpenoid).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a hopanoid, False otherwise.
        str: Explanation for the classification.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Get the Murcko scaffold, which represents the ring framework of the molecule.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    
    # Count the number of carbon atoms in the scaffold.
    scaffold_c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetSymbol() == "C")
    
    # Get the number of rings in the scaffold.
    ring_info = scaffold.GetRingInfo()
    num_rings = ring_info.NumRings() if ring_info is not None else 0
    
    # Count the total number of carbon atoms in the molecule as a sanity check.
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")
    
    # Check that the scaffold’s carbon count is in the expected range for a hopane skeleton.
    # Hopane itself is C30; allowing a slight margin here.
    if not (27 <= scaffold_c_count <= 33):
        return False, f"Murcko scaffold has {scaffold_c_count} carbons, not in the expected range for a hopane skeleton."
    
    # A hopane skeleton is pentacyclic (5 fused rings); check for at least 5 rings.
    if num_rings < 5:
        return False, f"Murcko scaffold has {num_rings} rings, fewer than 5 required for a hopane skeleton."
    
    # Finally, require that the molecule has enough carbons overall to be a triterpenoid.
    if total_c < 25:
        return False, "Total carbon count is too low for a triterpenoid."
    
    return True, "Molecule displays a fused polycyclic scaffold with ~30 carbons, consistent with a hopane-based triterpenoid (hopanoid)."

# Example usage:
if __name__ == "__main__":
    # You can test with one of the provided SMILES, for example dustanin:
    dustanin = "CC(C)(O)[C@H]1CC[C@@]2(C)[C@H]1C[C@H](O)[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3CC[C@@]12C"
    result, reason = is_hopanoid(dustanin)
    print("Is hopanoid:", result)
    print("Reason:", reason)