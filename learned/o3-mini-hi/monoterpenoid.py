"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: Monoterpenoid
Definition: Any terpenoid derived from a monoterpene. The term includes compounds in which the C10 skeleton of 
the parent monoterpene has been rearranged or modified by the removal of one or more skeletal atoms (generally methyl groups).

Heuristic improvements:
 1. Parse and sanitize the molecule.
 2. Obtain the Bemis–Murcko scaffold and break it into disconnected fragments.
 3. For each fragment, compute the number of carbon atoms, the total heavy atoms, and the ratio (carbon fraction).
 4. Also compute the aromatic fraction in the fragment – monoterpenoids are usually aliphatic.
 5. Prefer a fragment with high carbon fraction (≥0.60) and low aromatic fraction (<0.3). If none qualifies, take the fragment with the highest carbon count.
 6. If the chosen fragment’s carbon count is between 6 and 13 (inclusive), classify as monoterpenoid.
 7. Otherwise, return false with a reason indicating whether the scaffold carbon count is too low or too high.
 
Note: This is a heuristic approach that may misclassify borderline cases.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES.

    A monoterpenoid is derived from a monoterpene (the parent being a C10 skeleton, possibly rearranged or modified).
    This heuristic extracts the Bemis–Murcko scaffold fragments and evaluates each for:
      - Carbon atom count,
      - Carbon fraction (≥0.60), and
      - Low aromatic content (aromatic fraction < 0.3).
    We then classify as monoterpenoid if one of the fragments has a carbon count between 6 and 13 (inclusive).

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        (bool, str): Tuple of (True, reason) if classified as a monoterpenoid;
                     (False, reason) otherwise.
    """

    # Parse the SMILES; disable sanitization so we can catch errors manually.
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return False, "Invalid SMILES string."

    # Manually sanitize the molecule.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during sanitization: {e}"

    # Attempt kekulization; if it fails, log and continue.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass  # Not critical if kekulization fails.

    # Obtain the Bemis–Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error generating Murcko scaffold: {e}"

    if scaffold is None:
        return False, "Could not generate a scaffold."

    # Break the scaffold into disconnected fragments.
    frags = Chem.GetMolFrags(scaffold, asMols=True)
    if not frags:
        return False, "No scaffold fragments found."

    # Set our thresholds.
    desired_c_fraction = 0.60
    max_allowed_aromatic_fraction = 0.30
    best_frag = None
    best_carbon_count = 0
    best_c_fraction = 0.0

    # Evaluate each scaffold fragment.
    for frag in frags:
        atoms = frag.GetAtoms()
        heavy_atoms = [atom for atom in atoms if atom.GetAtomicNum() > 1]
        total_heavy = len(heavy_atoms)
        if total_heavy == 0:
            continue
        # Count carbons.
        c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        # Calculate carbon fraction.
        c_fraction = c_count / total_heavy

        # Compute aromatic atoms count.
        aromatic_atoms = [atom for atom in atoms if atom.GetIsAromatic()]
        aromatic_fraction = len(aromatic_atoms) / len(atoms) if atoms else 0

        # Prefer fragments that meet desired carbon fraction and low aromatic content.
        if c_fraction >= desired_c_fraction and aromatic_fraction < max_allowed_aromatic_fraction:
            # Look for the fragment with largest carbon count.
            if c_count > best_carbon_count:
                best_frag = frag
                best_carbon_count = c_count
                best_c_fraction = c_fraction
        # If none meet the criteria, we still keep track of the fragment with highest carbon count.
        elif best_frag is None and c_count > best_carbon_count:
            best_frag = frag
            best_carbon_count = c_count
            best_c_fraction = c_fraction

    if best_frag is None:
        return False, "No scaffold fragment selected."

    # Check the best fragment carbon count against the expected monoterpene range.
    # We allow a small relaxation at the lower end (6 instead of 7) to catch borderline cases.
    if best_carbon_count < 6:
        return False, f"Scaffold carbon count too low ({best_carbon_count}). Expected at least 6 carbons for a monoterpenoid core."
    if best_carbon_count > 13:
        return False, f"Scaffold carbon count too high ({best_carbon_count}). Expected 13 or fewer carbons for a monoterpenoid core."

    # Form a reason message.
    reason = (f"Scaffold carbon count is {best_carbon_count} with carbon fraction {best_c_fraction:.2f} "
              f"and low aromatic fraction (<{max_allowed_aromatic_fraction}); consistent with a monoterpenoid core.")
    return True, reason

# Example usage:
if __name__ == "__main__":
    # Test a few example SMILES
    examples = {
        "L-Menthone 1,2-glycerol ketal": "O1C2(OCC1CO)[C@@H](CC[C@H](C2)C)C(C)C",
        "oleuropein": "[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C",
        "xi-p-Mentha-3,8-dien-1-ol": "OC1(CCC(=CC1)C(C)=C)"
    }
    for name, smi in examples.items():
        classification, msg = is_monoterpenoid(smi)
        print(f"{name} --> {classification}: {msg}")