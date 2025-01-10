"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

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
    
    # Check for macrocyclic rings (12 or more members)
    ring_info = mol.GetRingInfo()
    macrocycle_found = False
    for atom_ring in ring_info.AtomRings():
        if len(atom_ring) >= 12:
            # Now, let's check if this ring has an ester linkage
            ring_atoms = set(atom_ring)
            found_ester = False
            for i, atom_idx in enumerate(atom_ring):
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 8:  # Check for Oxygen in the ring
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() in ring_atoms and neighbor.GetAtomicNum() == 6:  # Neighbor should be a Carbon
                            # Check if it's a carbonyl carbon (connected to oxygen through double bond)
                            for carbon_neighbor in neighbor.GetNeighbors():
                                if carbon_neighbor.GetAtomicNum() == 8 and neighbor.GetBondWith(atom).GetBondType() == Chem.BondType.SINGLE:
                                    if neighbor.GetBondWith(carbon_neighbor).GetBondType() == Chem.BondType.DOUBLE:
                                        found_ester = True
                                        break
                if found_ester:
                    macrocycle_found = True
                    break

            if macrocycle_found:
                break  # Break outer loop if macrocyclic lactone is found

    if not macrocycle_found:
        return False, "No macrocyclic lactone with 12 or more members found"

    # Assume additional polyketide structure features if macrocyclic lactone is present
    return True, "Contains a macrocyclic lactone with 12 or more members"