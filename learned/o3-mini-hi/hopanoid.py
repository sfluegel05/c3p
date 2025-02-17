"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: A hopanoid is defined as a triterpenoid based on a hopane skeleton.
This improved approach extracts the molecule's Murcko scaffold (selecting the largest 
fragment if disconnected), checks that:
  - the scaffold has at least 20 carbon atoms,
  - the scaffold has at least 5 rings,
  - the scaffold is mostly carbon (low heteroatom ratio),
and finally it requires that the scaffold contain a hopane‐like core (via a substructure match
to a hopane SMARTS pattern). These changes are designed to help reduce both the false negatives 
(from overly strict upper limits on carbon count) and the false positives (by ensuring the core 
is hopane‐like).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_hopanoid(smiles: str):
    """
    Determines whether a molecule is a hopanoid based on its SMILES string.
    A hopanoid is defined as a triterpenoid based on a hopane skeleton.
    This function:
      1. Extracts the molecule’s Murcko scaffold and (if needed) selects its largest fragment.
      2. Computes the number of carbon atoms and the heteroatom ratio in the scaffold.
      3. Checks that the scaffold has at least 5 rings.
      4. Requires that the overall molecule has at least 25 carbons.
      5. Checks that the scaffold matches a hopane-like core substructure.
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule is classified as a hopanoid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES to get an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Extract the Murcko scaffold.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Could not extract a Murcko scaffold."
    
    # If the scaffold is disconnected, select the largest fragment.
    frags = Chem.GetMolFrags(scaffold, asMols=True, sanitizeFrags=True)
    if frags:
        scaffold = max(frags, key=lambda m: m.GetNumAtoms())

    # Count carbons in the scaffold.
    scaffold_c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetSymbol() == "C")
    scaffold_total_atoms = scaffold.GetNumAtoms()
    non_c_count = scaffold_total_atoms - scaffold_c_count
    hetero_ratio = non_c_count / scaffold_total_atoms if scaffold_total_atoms > 0 else 1.0

    # Get ring information.
    ring_info = scaffold.GetRingInfo()
    num_rings = ring_info.NumRings() if ring_info is not None else 0

    # Total carbons in whole molecule.
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")

    # Heuristic checks:
    # (1) The scaffold must have at least 20 carbons (hopane cores are about C30 but modifications
    #     may add extra carbons; we use a minimum only).
    if scaffold_c_count < 20:
        return False, f"Murcko scaffold has {scaffold_c_count} carbons; expected at least 20 for a hopane core."
    
    # (2) At least 5 rings are required (the hopane skeleton is pentacyclic).
    if num_rings < 5:
        return False, f"Murcko scaffold has {num_rings} rings; at least 5 rings required."
    
    # (3) The scaffold should be almost all carbon (allowing a small margin, e.g. ≤15% heteroatoms).
    if hetero_ratio > 0.15:
        return False, (f"Murcko scaffold contains {non_c_count} non-carbon atoms out of {scaffold_total_atoms} "
                       f"({hetero_ratio*100:.1f}% heteroatoms); expected a primarily carbon scaffold.")
    
    # (4) Overall, the molecule should have enough carbons to be a triterpenoid.
    if total_c < 25:
        return False, "Total carbon count is too low for a triterpenoid."
    
    # (5) Check for a hopane-like core by substructure matching.
    # The following SMARTS string represents a simplified core of the hopane skeleton.
    hopane_template_smarts = "C1CC2CCC3C4CCC5CC4CC3CCC2C1"
    hopane_template = Chem.MolFromSmarts(hopane_template_smarts)
    if hopane_template is None:
        return False, "Internal error: invalid hopane template SMARTS."
    
    if not scaffold.HasSubstructMatch(hopane_template):
        return False, "Murcko scaffold does not contain a hopane-like core substructure."

    return True, ("Molecule has a fused polycyclic scaffold with a mostly carbon backbone "
                  f"({scaffold_c_count} C, {num_rings} rings) and a hopane core, consistent with a hopanoid.")

# Example usage:
if __name__ == "__main__":
    # Test with an example hopanoid (e.g., dustanin).
    dustanin = "CC(C)(O)[C@H]1CC[C@@]2(C)[C@H]1C[C@H](O)[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3CC[C@@]12C"
    result, reason = is_hopanoid(dustanin)
    print("Is hopanoid:", result)
    print("Reason:", reason)