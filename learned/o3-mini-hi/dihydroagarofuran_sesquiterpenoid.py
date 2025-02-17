"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: Dihydroagarofuran sesquiterpenoid
Definition: Any sesquiterpenoid with a dihydroagarofuran skeleton.
The dihydroagarofuran skeleton is expected to be a non‐aromatic (i.e. fully saturated) tricyclic core
with 15 carbon atoms (the standard sesquiterpene backbone) after removal of decorating substituents.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    The method extracts the Bemis–Murcko scaffold (i.e. the core of the molecule) and then checks:
      - That there are exactly 15 carbon atoms (sesquiterpenoid scaffold)
      - That the scaffold is tricyclic (the dihydroagarofuran core is composed of three fused rings)
      - That the scaffold contains no aromatic atoms (dihydro indicating saturation)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification (or non-classification)
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Extract the Bemis–Murcko scaffold (removes terminal substituents)
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Could not extract scaffold from molecule."

    # Count carbon atoms in the scaffold.
    carbon_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Get number of rings in the scaffold.
    ring_info = scaffold.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    # Check for aromatic atoms in the scaffold.
    aromatic_atoms = [atom for atom in scaffold.GetAtoms() if atom.GetIsAromatic()]

    # For a dihydroagarofuran skeleton we expect:
    # - Exactly 15 carbons in the core backbone (sesquiterpenoid backbone)
    # - The core is built from three fused rings (tricyclic system)
    # - The core is fully saturated (no aromatic atoms present)
    if carbon_count != 15:
        return (False, f"Scaffold has {carbon_count} carbon atoms; expected 15 for a sesquiterpenoid skeleton.")
    if num_rings != 3:
        return (False, f"Scaffold has {num_rings} rings; expected a tricyclic (3-ring) system for dihydroagarofuran.")
    if aromatic_atoms:
        return (False, "Scaffold contains aromatic atoms; expected a fully saturated (dihydro) skeleton.")
    
    return (True, "Molecule contains a 15-carbon, non-aromatic, tricyclic scaffold consistent with dihydroagarofuran sesquiterpenoid.")

# Example usage (uncomment to test):
# smiles = "[C@@]1(C(O[C@@]2([C@](O)([C@@]34[C@H](OC(C)=O)..."
# result, reason = is_dihydroagarofuran_sesquiterpenoid(smiles)
# print(result, reason)