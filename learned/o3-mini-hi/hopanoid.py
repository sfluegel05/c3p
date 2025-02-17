"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: A hopanoid is defined as a triterpenoid based on a hopane skeleton.
This implementation processes the molecule by extracting its Murcko scaffold (selecting the largest 
fragment if disconnected), checking that:
  - the scaffold has at least 20 carbon atoms,
  - the scaffold has at least 5 rings,
  - the scaffold is almost entirely carbon (≤15% heteroatoms),
and then requires that the scaffold match a hopane‐like core (via a substructure match to a simplified
hopane SMARTS pattern). Additionally, all primary RDKit calls are wrapped in try/except blocks to catch
kekulization or sanitization issues.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_hopanoid(smiles: str):
    """
    Determines whether a molecule is a hopanoid based on its SMILES string.
    A hopanoid is defined as a triterpenoid based on a hopane skeleton.
    This function performs the following:
      1. Parses the SMILES string (with error handling) and ensures proper kekulization.
      2. Extracts the Murcko scaffold; if the scaffold is made of disconnected fragments, selects the largest.
      3. Evaluates the scaffold for the following:
            - At least 20 carbons.
            - At least 5 rings (hopane is pentacyclic).
            - A low heteroatom ratio (≤15% non-carbon atoms).
      4. Ensures that the full molecule has at least 25 carbon atoms.
      5. Searches for a simplified hopane-like core within the scaffold using a SMARTS substructure.
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule qualifies as a hopanoid, False otherwise.
        str: Explanation for the classification.
    """
    try:
        # Parse the SMILES into an RDKit molecule (includes sanitization)
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception as e:
        return False, f"Error parsing SMILES: {str(e)}"
    
    if mol is None:
        return False, "Invalid SMILES string."

    # Attempt to kekulize the molecule (this can fail for problematic molecules)
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        return False, f"Molecule could not be kekulized: {str(e)}"

    # Extract the Murcko scaffold (wrap in try/except in case of internal errors)
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error extracting Murcko scaffold: {str(e)}"
    
    if scaffold is None:
        return False, "Could not extract a Murcko scaffold."
    
    # If the scaffold is disconnected, select the largest fragment.
    frags = Chem.GetMolFrags(scaffold, asMols=True, sanitizeFrags=True)
    if frags:
        scaffold = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Count the number of carbon atoms in the scaffold.
    scaffold_c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetSymbol() == "C")
    scaffold_total_atoms = scaffold.GetNumAtoms()
    non_c_count = scaffold_total_atoms - scaffold_c_count
    hetero_ratio = non_c_count / scaffold_total_atoms if scaffold_total_atoms > 0 else 1.0
    
    # Count the number of rings in the scaffold.
    ring_info = scaffold.GetRingInfo()
    num_rings = ring_info.NumRings() if ring_info is not None else 0

    # Total carbon count in the whole molecule.
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")
    
    # Apply heuristic thresholds:
    if scaffold_c_count < 20:
        return False, f"Murcko scaffold has {scaffold_c_count} carbons; expected at least 20 for a hopane core."
    if num_rings < 5:
        return False, f"Murcko scaffold has {num_rings} rings; at least 5 rings required (hopane is pentacyclic)."
    if hetero_ratio > 0.15:
        return False, (f"Murcko scaffold contains {non_c_count} non-carbon atoms out of {scaffold_total_atoms} "
                       f"({hetero_ratio*100:.1f}% heteroatoms); expected a predominantly carbon scaffold.")
    if total_c < 25:
        return False, "Total carbon count is too low for a triterpenoid."
    
    # Create a simplified hopane-like core template via SMARTS.
    hopane_template_smarts = "C1CC2CCC3C4CCC5CC4CC3CCC2C1"
    hopane_template = Chem.MolFromSmarts(hopane_template_smarts)
    if hopane_template is None:
        return False, "Internal error: invalid hopane template SMARTS."
    
    if not scaffold.HasSubstructMatch(hopane_template):
        return False, "Murcko scaffold does not contain a hopane-like core substructure."
    
    return True, (f"Molecule qualifies as a hopanoid with a scaffold featuring {scaffold_c_count} carbons "
                  f"and {num_rings} rings, consistent with a hopane skeleton.")

# Example usage:
if __name__ == "__main__":
    # Test with a reported hopanoid (e.g., dustanin).
    dustanin = "CC(C)(O)[C@H]1CC[C@@]2(C)[C@H]1C[C@H](O)[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3CC[C@@]12C"
    result, explanation = is_hopanoid(dustanin)
    print("Is hopanoid:", result)
    print("Explanation:", explanation)