"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: Monoterpenoid (any terpenoid derived from a monoterpene, with a C10-derived skeleton possibly rearranged or modified)

Heuristic:
 1. Parse the molecule and sanitize.
 2. Attempt to kekulize the molecule (clearing aromatic flags) to avoid kekulization errors.
 3. Obtain the Bemis–Murcko scaffold of the molecule.
 4. Break the scaffold into disconnected fragments.
 5. For each fragment, compute the number of carbons and the ratio of carbon atoms to heavy atoms.
 6. Select the fragment with a high carbon fraction (>=0.60) and the highest carbon count. If none meet the threshold, choose the one with the highest carbon count.
 7. If the chosen fragment has between 7 and 13 carbons (inclusive) the molecule is classified as a monoterpenoid.
 
If any error occurs (e.g., during kekulization or scaffold generation), we now catch the exception and return a helpful message.
 
Examples of structures that belong to this class are provided in the task description.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    
    This function uses the Bemis–Murcko scaffold to isolate the core of the molecule,
    breaks it into fragments and then measures the carbon content. A candidate fragment
    is selected based on its carbon-to-heavy-atom ratio (>= 0.60) and its total number 
    of carbon atoms (expected to be between 7 and 13 for a monoterpenoid).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if the molecule is classified 
                     as a monoterpenoid, False otherwise, along with a reason.
    """
    # Parse the SMILES into an RDKit molecule. We disable automatic sanitization to handle kekulization manually.
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return False, "Invalid SMILES string."

    # Attempt manual sanitization.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during sanitization: {e}"

    # Try to kekulize with aromatic flags cleared, to avoid errors with unusual valences.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        # Instead of failing, we log the error and continue with the molecule as is.
        # This approach can help process molecules with non-standard aromaticity.
        pass

    # Try to obtain the Bemis–Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error generating Murcko scaffold: {e}"

    if scaffold is None:
        return False, "Could not generate a scaffold."

    # Split the scaffold into disconnected fragments.
    frags = Chem.GetMolFrags(scaffold, asMols=True)
    if not frags:
        return False, "No scaffold fragments found."

    desired_ratio = 0.60
    best_frag = None
    best_carbon_count = 0
    best_ratio = 0.0

    # Process each fragment to compute carbon counts and ratios.
    for frag in frags:
        atoms = frag.GetAtoms()
        total_heavy = sum(1 for atom in atoms if atom.GetAtomicNum() > 1)  # only heavy atoms
        carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if total_heavy == 0:
            continue
        ratio = carbon_count / total_heavy
        # Select fragment if it meets the desired ratio threshold and has more carbons than previous.
        if ratio >= desired_ratio and carbon_count > best_carbon_count:
            best_frag = frag
            best_carbon_count = carbon_count
            best_ratio = ratio

    # If no fragment met the ratio requirement, choose the one with the highest carbon count overall.
    if best_frag is None:
        for frag in frags:
            atoms = frag.GetAtoms()
            carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
            if carbon_count > best_carbon_count:
                best_frag = frag
                best_carbon_count = carbon_count
        best_ratio = None  # did not meet the desired ratio

    # Final check based on the number of carbons in the selected fragment.
    if best_carbon_count < 7:
        return False, f"Scaffold carbon count too low ({best_carbon_count}). Expected at least 7 carbons."
    elif best_carbon_count > 13:
        return False, f"Scaffold carbon count too high ({best_carbon_count}). Expected 13 or fewer carbons."
    
    # Return successful classification with reasoning.
    if best_ratio is not None:
        reason = (f"Scaffold carbon count is {best_carbon_count} with high carbon fraction (≥{desired_ratio}); "
                  "consistent with a monoterpenoid core.")
    else:
        reason = (f"Scaffold carbon count is {best_carbon_count} (selected fragment did not meet desired carbon fraction threshold); "
                  "borderline monoterpenoid core.")
    
    return True, reason