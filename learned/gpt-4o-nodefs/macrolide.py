"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    Macrolides are characterized by a large macrocyclic lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Identify candidate macrolide rings
    for ring in atom_rings:
        # Check if ring size is in typical macrolide range
        if 12 <= len(ring) <= 16:
            # Check for an ester linkage within the ring
            found_ester = False
            for bond in ring:
                atom1 = mol.GetAtomWithIdx(bond)
                for neighbor in atom1.GetNeighbors():
                    # Identify the ester pattern: "C(=O)O"
                    if neighbor.GetAtomicNum() == 8:  # Oxygen
                        carbon_atom = neighbor.GetNeighbors()[0]
                        if carbon_atom.GetAtomicNum() == 6:  # Check carbon
                            if carbon_atom.GetTotalNumHs() == 0 and any(nb.GetAtomicNum() == 8 for nb in carbon_atom.GetNeighbors() if nb.GetIdx() != atom1.GetIdx()):
                                found_ester = True
                                break
                if found_ester:
                    break
            if found_ester:
                return True, "Contains macrocyclic lactone ring with an ester linkage"

    return False, "No characteristic macrolide macrocyclic lactone ring found"