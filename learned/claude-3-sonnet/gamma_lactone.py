"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: CHEBI:35524 gamma-lactone
A lactone having a five-membered lactone ring.
"""
from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined as a five-membered lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    
    # Iterate over rings
    for ring in ring_info.AtomRings():
        if len(ring) == 5:  # Five-membered ring
            oxygen_atoms = [mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in ring]
            if sum(oxygen_atoms) == 1:  # Exactly one oxygen atom
                oxygen_idx = ring[oxygen_atoms.index(True)]
                oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
                
                # Check if oxygen is part of a lactone
                if oxygen_atom.IsInRingSize(5) and oxygen_atom.GetIsAromatic() == False:
                    neighbors = [mol.GetAtomWithIdx(neighbor_idx) for neighbor_idx in oxygen_atom.GetNeighbors()]
                    carbonyl_carbons = [neighbor for neighbor in neighbors if neighbor.GetSymbol() == 'C' and neighbor.GetFormalCharge() == 0]
                    if len(carbonyl_carbons) == 1:
                        carbonyl_carbon = carbonyl_carbons[0]
                        neighbors_of_carbonyl = [mol.GetAtomWithIdx(neighbor_idx) for neighbor_idx in carbonyl_carbon.GetNeighbors()]
                        if len([neighbor for neighbor in neighbors_of_carbonyl if neighbor.GetSymbol() == 'O' and neighbor.GetFormalCharge() == 0]) == 1:
                            return True, "Contains a five-membered lactone ring"
    
    return False, "No gamma-lactone substructure found"