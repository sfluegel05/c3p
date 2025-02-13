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
    is attached to a carbon atom that is itself bonded to an aromatic ring.

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
    aromatic_rings = mol.GetAromaticRings()
    
    # Check if any carbon with a primary alcoholic group is attached to an aromatic ring
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetTotalNumHs() == 2:  # primary carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == "O" and neighbor.GetDegree() == 1:  # hydroxyl oxygen
                    # Check if the primary carbon is attached to an aromatic ring
                    for ring in aromatic_rings:
                        if atom.IsMemberOfAromaticRing(ring):
                            return True, "Molecule contains an aromatic primary alcohol group"

    return False, "Molecule does not contain an aromatic primary alcohol group"