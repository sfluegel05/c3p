"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: CHEBI:35476 aromatic primary alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol is any primary alcohol where the alcoholic hydroxy group
    is directly attached to a carbon atom that is part of an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find aromatic rings
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if mol.RingIsPlanar(ring)]
    
    # Check if any aromatic ring has a carbon with a primary alcoholic group attached
    for ring in aromatic_rings:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetHybridization() == Chem.HybridizationType.SP2 and atom.GetTotalNumHs() == 0:  # aromatic carbon
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == "O" and neighbor.GetDegree() == 1:  # hydroxyl oxygen
                        for carbon_neighbor in neighbor.GetNeighbors():
                            if carbon_neighbor.GetSymbol() == "C" and carbon_neighbor.GetTotalNumHs() == 2:  # primary carbon
                                return True, "Molecule contains an aromatic primary alcohol group"

    return False, "Molecule does not contain an aromatic primary alcohol group"