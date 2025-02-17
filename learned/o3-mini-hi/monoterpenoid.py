"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: Monoterpenoid

Definition: Any terpenoid derived from a monoterpene. The term includes compounds in which
the C10 skeleton of the parent monoterpene has been rearranged or modified (often with some
methyl groups removed). As monoterpenoids tend to be aliphatic with a small core (roughly 6–13 carbons)
we use a heuristic based on the Bemis–Murcko scaffold fragments. For each fragment we compute:
  – Carbon count and total heavy atoms (from which a carbon fraction is obtained),
  – Aromatic atom fraction, and
  – An oxygen-to-carbon ratio (to help screen out sugar‐like fragments).
A fragment with high carbon fraction (≥0.75), low aromatic fraction (<0.3) and O/C ratio (<0.8)
is selected. If its carbon count falls between 6 and 13, the molecule is classified as a monoterpenoid.
This heuristic may miss borderline cases.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES.

    The approach is as follows:
      1. Parse and sanitize the molecule.
      2. Obtain the Bemis–Murcko scaffold and break it into disconnected fragments.
      3. For each fragment, compute:
           - The number of carbon atoms,
           - The total number of heavy (non-hydrogen) atoms,
           - Carbon fraction (carbon count/total heavy atoms),
           - Aromatic fraction, and
           - Oxygen-to-carbon ratio.
      4. Select fragments that have a high carbon fraction (>=0.75),
         low aromatic fraction (<0.3), and low oxygen-to-carbon ratio (<0.8).
      5. If at least one fragment has a carbon count between 6 and 13 (inclusive),
         classify as monoterpenoid.
      6. Otherwise return False with a reason indicating whether the selected fragment
         has too few or too many carbons, or if no suitable scaffold fragment was found.

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns: 
        (bool, str): Tuple of (True, reason) if classified as a monoterpenoid;
                     otherwise (False, reason).
    """
    # Parse the SMILES with sanitization disabled (we will do it manually)
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return False, "Invalid SMILES string."

    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during sanitization: {e}"

    # Try kekulization; if it fails, we continue.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass

    # Generate the Bemis–Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error generating Murcko scaffold: {e}"

    if scaffold is None:
        return False, "Could not generate a scaffold."

    # Break the scaffold into disconnected fragments.
    fragments = Chem.GetMolFrags(scaffold, asMols=True)
    if not fragments:
        return False, "No scaffold fragments found."

    # Set thresholds.
    desired_c_fraction = 0.75
    max_allowed_arom_fraction = 0.30
    max_allowed_o_to_c = 0.80

    best_frag = None
    best_carbon_count = 0
    best_c_fraction = 0.0

    # Evaluate each scaffold fragment.
    for frag in fragments:
        atoms = list(frag.GetAtoms())
        # Select heavy (non-hydrogen) atoms.
        heavy_atoms = [atom for atom in atoms if atom.GetAtomicNum() > 1]
        if not heavy_atoms:
            continue
        total_heavy = len(heavy_atoms)
        # Count carbon atoms.
        c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        # Count oxygen atoms.
        o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        c_fraction = c_count / total_heavy

        # Count aromatic atoms.
        aromatic_atoms = [atom for atom in atoms if atom.GetIsAromatic()]
        arom_fraction = len(aromatic_atoms) / len(atoms) if atoms else 0

        # Calculate oxygen-to-carbon ratio (avoid division by zero).
        o_to_c_ratio = o_count / c_count if c_count > 0 else 0

        # Prefer fragments that meet desired thresholds:
        if c_fraction >= desired_c_fraction and arom_fraction < max_allowed_arom_fraction and o_to_c_ratio < max_allowed_o_to_c:
            # Choose the fragment with the highest carbon count that satisfies our criteria.
            if c_count > best_carbon_count:
                best_frag = frag
                best_carbon_count = c_count
                best_c_fraction = c_fraction
        # If none meet the strict criteria so far, keep track of a fallback fragment.
        elif best_frag is None and c_count > best_carbon_count:
            best_frag = frag
            best_carbon_count = c_count
            best_c_fraction = c_fraction

    if best_frag is None:
        return False, "No scaffold fragment selected."

    # Check if the selected fragment carbon count lies within expected bounds for a monoterpenoid core.
    if best_carbon_count < 6:
        return False, f"Scaffold carbon count too low ({best_carbon_count}). Expected at least 6 carbons for a monoterpenoid core."
    if best_carbon_count > 13:
        return False, f"Scaffold carbon count too high ({best_carbon_count}). Expected 13 or fewer carbons for a monoterpenoid core."

    reason = (f"Scaffold carbon count is {best_carbon_count} with carbon fraction {best_c_fraction:.2f} "
              f"and low aromatic fraction (<{max_allowed_arom_fraction}); consistent with a monoterpenoid core.")
    return True, reason

# Example usage:
if __name__ == "__main__":
    examples = {
        "L-Menthone 1,2-glycerol ketal": "O1C2(OCC1CO)[C@@H](CC[C@H](C2)C)C(C)C",
        "oleuropein": "[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C",
        "xi-p-Mentha-3,8-dien-1-ol": "OC1(CCC(=CC1)C(C)=C)"
    }
    for name, smi in examples.items():
        classification, msg = is_monoterpenoid(smi)
        print(f"{name}: {classification}\nReason: {msg}\n")