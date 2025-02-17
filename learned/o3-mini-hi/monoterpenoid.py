"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: Monoterpenoid (any terpenoid derived from a monoterpene, with a C10-derived skeleton that may be rearranged or modified)

Heuristic improvements:
  1. Parse the molecule and get its Bemis–Murcko scaffold.
  2. Break the scaffold into disconnected fragments.
  3. For each fragment compute two values:
       - The number of carbon atoms.
       - The ratio of carbon atoms to total heavy atoms (to discount fragments such as sugars with many oxygens).
  4. Select the fragment that has a high carbon fraction (we require ≥ 0.60) and has the highest carbon count.
  5. If no fragment meets the desired ratio, fall back on the fragment with the highest carbon count.
  6. Finally, if the chosen fragment has between 7 and 13 carbons (inclusive) we classify the molecule as a monoterpenoid;
     otherwise it is rejected.
Note: This is a heuristic classification and may result in false positives or negatives.
     
Examples: See provided outcomes.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    
    The method uses the Bemis–Murcko scaffold of the molecule and then breaks it into fragments.
    It selects the fragment with a sufficiently high carbon-to-heavy-atom ratio (>=0.60) and the 
    largest number of carbons. If the number of carbons in this candidate fragment is between 7 and 13,
    it is considered consistent with a monoterpenoid core.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element indicates if the molecule is classified as a monoterpenoid,
                     and the second element is the reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    try:
        # Obtain the Bemis–Murcko scaffold.
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error generating Murcko scaffold: {e}"
    
    if scaffold is None:
        return False, "Could not generate a scaffold."
    
    # Break the scaffold into disconnected fragments.
    frags = Chem.GetMolFrags(scaffold, asMols=True)
    if not frags:
        return False, "No scaffold fragments found."
    
    # Define threshold for the ratio of carbons to heavy atoms.
    desired_ratio = 0.60
    best_frag = None
    best_carbon_count = 0
    best_ratio = 0.0
    
    # Process each fragment.
    for frag in frags:
        atoms = frag.GetAtoms()
        total_heavy = sum(1 for atom in atoms if atom.GetAtomicNum() > 1)  # ignoring hydrogens
        carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if total_heavy == 0:
            continue
        ratio = carbon_count / total_heavy
        # If the fragment has a high carbon fraction and more carbons than previously seen, choose it.
        if ratio >= desired_ratio and carbon_count > best_carbon_count:
            best_frag = frag
            best_carbon_count = carbon_count
            best_ratio = ratio
    # If none of the fragments pass the ratio threshold, then simply pick the fragment with the highest carbon count.
    if best_frag is None:
        for frag in frags:
            atoms = frag.GetAtoms()
            total_heavy = sum(1 for atom in atoms if atom.GetAtomicNum() > 1)
            carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
            if carbon_count > best_carbon_count:
                best_frag = frag
                best_carbon_count = carbon_count
        best_ratio = None  # ratio check was not met
    
    # Final classification based on the chosen fragment's carbon count.
    if best_carbon_count < 7:
        return False, f"Scaffold carbon count too low ({best_carbon_count}). Expected at least 7 carbons."
    elif best_carbon_count > 13:
        return False, f"Scaffold carbon count too high ({best_carbon_count}). Expected 13 or fewer carbons."
    
    # Build reason message.
    if best_ratio is not None:
        reason = f"Scaffold carbon count is {best_carbon_count} with high carbon fraction (≥{desired_ratio}); consistent with a monoterpenoid core."
    else:
        reason = f"Scaffold carbon count is {best_carbon_count} (selected fragment did not meet desired carbon fraction threshold); borderline monoterpenoid core."
    
    return True, reason