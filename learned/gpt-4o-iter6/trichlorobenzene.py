"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is a member of the chlorobenzenes class with three chloro substituents on at least one benzene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find aromatic rings in the molecule
    aromatic_rings = mol.GetAromaticRings()
    
    # Check each aromatic ring for Cl substituents
    for ring in aromatic_rings:
        chloro_count = sum(1 for atom_idx in ring 
                           if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C' 
                           and any(neighbor.GetSymbol() == 'Cl' for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors()))
        
        if chloro_count == 3:
            return True, "Contains a benzene ring with three chloro substituents"
    
    return False, "Benzene ring with three chloro substituents not found"