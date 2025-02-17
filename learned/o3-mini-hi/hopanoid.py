"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: A hopanoid is defined as a triterpenoid based on a hopane skeleton.
This implementation:
  1. Parses the SMILES and kekulizes the molecule.
  2. Extracts the Murcko scaffold (using the largest fragment if disconnected).
  3. Checks that the scaffold has at least 20 carbon atoms, at least 5 rings, and is nearly all carbon (≤15% heteroatoms).
  4. Checks that the entire molecule has at least 25 carbons.
  5. Finally, requires that the scaffold contain a hopane‐like core via a SMARTS substructure match.
Note: The chosen hopane template SMARTS is a simplified pattern for a fused pentacyclic ring system.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_hopanoid(smiles: str):
    """
    Determines whether a molecule is a hopanoid based on its SMILES string.
    A hopanoid is defined as a triterpenoid based on a hopane skeleton.
    
    The function:
      1. Parses the SMILES (with error handling and kekulization).
      2. Extracts the Murcko scaffold; if disconnected, selects the largest fragment.
      3. Ensures that the scaffold has at least 20 carbons, at least 5 rings, and a low heteroatom ratio.
      4. Checks the full molecule contains at least 25 carbon atoms.
      5. Uses a simplified hopane SMARTS pattern (a fused pentacyclic core) to search within the scaffold.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple where the first item is True if the molecule qualifies as a hopanoid,
                     and the second is an explanation.
                     On error or inability to classify, returns (False, <explanation>).
    """
    # Step 1. Parse the SMILES string with sanitization.
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception as e:
        return False, f"Error parsing SMILES: {str(e)}"
    if mol is None:
        return False, "Invalid SMILES string."

    # Step 2. Kekulize the molecule (kekulization may fail for some structures).
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        return False, f"Molecule could not be kekulized: {str(e)}"

    # Step 3. Extract the Murcko scaffold (may raise an error if problematic).
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

    # Step 4. Count carbon atoms and rings in the scaffold.
    scaffold_c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetSymbol() == "C")
    total_scaffold_atoms = scaffold.GetNumAtoms()
    non_c_count = total_scaffold_atoms - scaffold_c_count
    hetero_ratio = non_c_count / total_scaffold_atoms if total_scaffold_atoms > 0 else 1.0

    ring_info = scaffold.GetRingInfo()
    num_rings = ring_info.NumRings() if ring_info is not None else 0

    # Also count carbons in the full molecule.
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")

    # Step 5. Apply heuristic thresholds.
    if scaffold_c_count < 20:
        return False, f"Murcko scaffold has {scaffold_c_count} carbons; expected at least 20 for a hopane core."
    if num_rings < 5:
        return False, f"Murcko scaffold has {num_rings} rings; at least 5 rings required (hopane is pentacyclic)."
    if hetero_ratio > 0.15:
        return False, (f"Murcko scaffold contains {non_c_count} non-carbon atoms out of {total_scaffold_atoms} "
                       f"({hetero_ratio*100:.1f}% heteroatoms); expected a predominantly carbon scaffold.")
    if total_c < 25:
        return False, "Total carbon count in the molecule is too low for a triterpenoid."

    # Step 6. Create a simplified hopane-like SMARTS pattern.
    # This template is designed to capture a fused pentacyclic system.
    hopane_template_smarts = "C1CCC2C3CCC4C5CCC4C3CC2C1"
    try:
        hopane_template = Chem.MolFromSmarts(hopane_template_smarts)
    except Exception as e:
        return False, f"Error in hopane template SMARTS: {str(e)}"
    if hopane_template is None:
        return False, "Internal error: invalid hopane template SMARTS."

    # Check if the scaffold contains the hopane-like core.
    if not scaffold.HasSubstructMatch(hopane_template):
        return False, "Murcko scaffold does not contain a hopane-like (fused pentacyclic) core substructure."

    # If all checks pass, classify the molecule as a hopanoid.
    return True, (f"Molecule qualifies as a hopanoid: scaffold contains {scaffold_c_count} carbons, "
                  f"{num_rings} rings, and a fused pentacyclic core consistent with a hopane skeleton.")

# Example usage (for testing):
if __name__ == "__main__":
    # Example: dustanin (a reported hopanoid)
    dustanin = "CC(C)(O)[C@H]1CC[C@@]2(C)[C@H]1C[C@H](O)[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3CC[C@@]12C"
    result, explanation = is_hopanoid(dustanin)
    print("Is hopanoid:", result)
    print("Explanation:", explanation)