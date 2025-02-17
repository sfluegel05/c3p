"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: A hopanoid is defined as a triterpenoid based on a hopane skeleton.
This implementation:
  1. Parses and kekulizes the molecule.
  2. Extracts the Murcko scaffold (choosing the largest fragment if disconnected).
  3. Checks that the scaffold has at least 20 carbon atoms, at least 5 rings, and a low heteroatom ratio.
  4. Checks that the full molecule has at least 25 carbons.
  5. Instead of a fixed SMARTS pattern (which proved too rigid), we require that the scaffold
     contain at least 5 six‐membered rings. This is in line with the hopane (pentacyclic) core.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_hopanoid(smiles: str):
    """
    Determines whether a molecule is a hopanoid (a triterpenoid based on a hopane skeleton)
    based on its SMILES string.
    
    The process is:
      1. Parse the SMILES (with sanitization and kekulization).
      2. Extract the Murcko scaffold (if disconnected, uses the largest fragment).
      3. Verify that the scaffold has at least 20 carbon atoms, at least 5 rings, and a low heteroatom ratio.
      4. Ensure that the full molecule has at least 25 carbon atoms.
      5. Analyze the Murcko scaffold’s ring system:
         - Count the rings that are exactly six-membered (cyclohexane type);
         - Require that at least 5 such rings exist to reflect the fused pentacyclic hopane core.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where the first element is True if the molecule qualifies as a hopanoid,
                     otherwise False. The second element explains the reasoning.
    """
    # Step 1: Parse the SMILES string.
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception as e:
        return False, f"Error parsing SMILES: {str(e)}"
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Step 2: Kekulize the molecule (this may fail for some structures).
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        return False, f"Molecule could not be kekulized: {str(e)}"
    
    # Step 3: Extract the Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error extracting Murcko scaffold: {str(e)}"
    if scaffold is None:
        return False, "Could not extract a Murcko scaffold."
    
    # If scaffold is disconnected, select the largest fragment.
    frags = Chem.GetMolFrags(scaffold, asMols=True, sanitizeFrags=True)
    if frags:
        scaffold = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Step 4: Verify scaffold carbon content and ring number.
    scaffold_c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetSymbol() == "C")
    total_scaffold_atoms = scaffold.GetNumAtoms()
    non_c_count = total_scaffold_atoms - scaffold_c_count
    hetero_ratio = non_c_count / total_scaffold_atoms if total_scaffold_atoms > 0 else 1.0

    # Get ring information for the scaffold.
    ring_info = scaffold.GetRingInfo()
    num_rings = ring_info.NumRings() if ring_info is not None else 0
    
    # Count carbons in the full molecule.
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")
    
    # Heuristic thresholds.
    if scaffold_c_count < 20:
        return False, f"Murcko scaffold has only {scaffold_c_count} carbons; expected at least 20 for a hopane core."
    if num_rings < 5:
        return False, f"Murcko scaffold has {num_rings} rings; at least 5 rings (pentacyclic) are required."
    if hetero_ratio > 0.15:
        return False, (f"Murcko scaffold has {non_c_count} non-carbon atoms out of {total_scaffold_atoms} "
                       f"({hetero_ratio*100:.1f}% heteroatoms); expected a predominantly carbon scaffold.")
    if total_c < 25:
        return False, "Total carbon count in the molecule is too low for a triterpenoid hopanoid."
    
    # Step 5: Analyze the ring system of the scaffold.
    # We require that the fused core contains at least 5 six-membered rings.
    rings = ring_info.AtomRings()
    six_membered_ring_count = sum(1 for r in rings if len(r) == 6)
    if six_membered_ring_count < 5:
        return False, (f"Scaffold contains only {six_membered_ring_count} six-membered rings; "
                       "expected at least 5 for a hopane (pentacyclic) core.")
    
    # If all heuristics pass, classify as hopanoid.
    return True, (f"Molecule qualifies as a hopanoid: scaffold has {scaffold_c_count} carbons, "
                  f"{num_rings} rings (including {six_membered_ring_count} six-membered rings), "
                  f"and overall {total_c} carbons in the molecule.")

# Example usage:
if __name__ == "__main__":
    # Test example: dustanin (a reported hopanoid)
    dustanin_smiles = "CC(C)(O)[C@H]1CC[C@@]2(C)[C@H]1C[C@H](O)[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3CC[C@@]12C"
    result, explanation = is_hopanoid(dustanin_smiles)
    print("Input SMILES:", dustanin_smiles)
    print("Is hopanoid:", result)
    print("Explanation:", explanation)