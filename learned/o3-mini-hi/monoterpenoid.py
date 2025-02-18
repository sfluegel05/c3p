"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: Monoterpenoid (any terpenoid derived from a monoterpene, with a C10‐derived skeleton possibly rearranged or modified)

Heuristic improvement:
 1. Parse and sanitize the molecule.
 2. Check the overall carbon fraction – many monoterpenoids (being non‐polar) have a high proportion of carbons.
 3. Obtain the Bemis–Murcko scaffold and break it into disconnected fragments.
 4. For each fragment, compute the number of carbon atoms, the number of heavy atoms, and their ratio.
 5. Prefer a fragment with a high carbon-to-heavy-atom ratio (≥0.60). If none qualify, take the fragment with the largest carbon count.
 6. If the best fragment contains between 7 and 13 carbons (inclusive), classify as monoterpenoid.
 7. Otherwise, return false with a reason stating whether the scaffold is too low or too high in carbon content.
 
Note: This is a heuristic approach and may misclassify borderline cases.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES.
    
    A monoterpenoid is derived from a monoterpene (i.e. a C10 skeleton, possibly rearranged or slightly modified)
    and is expected to have a Bemis–Murcko scaffold with a high carbon fraction. In this heuristic, we
    require that (A) the overall molecule is carbon-rich (≥75% of its heavy atoms are carbon), and (B)
    one of the scaffold fragments has between 7 and 13 carbon atoms and a high carbon fraction (≥60%).
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        (bool, str): Tuple of (True, reason) if classified as a monoterpenoid;
                     (False, reason) otherwise.
    """
    
    # Parse the SMILES (disable automatic sanitization to perform manual checks)
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return False, "Invalid SMILES string."
    
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during sanitization: {e}"
    
    # Attempt kekulization; if it fails, log and continue.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass  # Continue even if kekulization fails
    
    # Additional overall check: the molecule (as a whole) should be carbon-rich.
    atoms = mol.GetAtoms()
    heavy_atoms = [atom for atom in atoms if atom.GetAtomicNum() > 1]
    total_heavy = len(heavy_atoms)
    total_carbons = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    if total_heavy == 0:
        return False, "No heavy atoms found in molecule."
    overall_c_ratio = total_carbons / total_heavy
    if overall_c_ratio < 0.75:
        return False, (f"Overall molecular carbon fraction too low ({overall_c_ratio:.2f}); "
                       "likely not a carbon‐rich monoterpenoid derivative.")
    
    # Obtain the Bemis–Murcko scaffold
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
    
    desired_ratio = 0.60
    best_frag = None
    best_carbon_count = 0
    best_ratio = 0.0
    
    # Evaluate each scaffold fragment.
    for frag in frags:
        frag_atoms = frag.GetAtoms()
        frag_heavy = [atom for atom in frag_atoms if atom.GetAtomicNum() > 1]
        frag_heavy_count = len(frag_heavy)
        frag_carbons = sum(1 for atom in frag_atoms if atom.GetAtomicNum() == 6)
        if frag_heavy_count == 0:
            continue
        frag_ratio = frag_carbons / frag_heavy_count
        
        # Prefer fragments meeting our desired carbon fraction threshold.
        if frag_ratio >= desired_ratio and frag_carbons > best_carbon_count:
            best_frag = frag
            best_carbon_count = frag_carbons
            best_ratio = frag_ratio
    
    # If no fragment meets the desired ratio threshold, simply take the fragment with the highest carbon count.
    if best_frag is None:
        for frag in frags:
            frag_carbons = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
            if frag_carbons > best_carbon_count:
                best_frag = frag
                best_carbon_count = frag_carbons
        best_ratio = None  # did not meet desired ratio threshold
    
    # Heuristic classification based on the chosen fragment's carbon count.
    if best_carbon_count < 7:
        return False, f"Scaffold carbon count too low ({best_carbon_count}). Expected at least 7 carbons."
    elif best_carbon_count > 13:
        return False, f"Scaffold carbon count too high ({best_carbon_count}). Expected 13 or fewer carbons."
    
    # If we reach here, the scaffold appears appropriate.
    if best_ratio is not None:
        reason = (f"Scaffold carbon count is {best_carbon_count} with high carbon fraction (≥{desired_ratio}); "
                  "consistent with a monoterpenoid core.")
    else:
        reason = (f"Scaffold carbon count is {best_carbon_count} (selected fragment did not meet desired carbon fraction threshold); "
                  "borderline monoterpenoid core.")
    
    return True, reason