"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    Macrolides are characterized by a large macrocyclic lactone ring, typically with an ester linkage.

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
    
    # Get the molecular rings information with a focus on larger macrocycles
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # SMARTS pattern for a lactone, which is a cyclic ester
    lactone_pattern = Chem.MolFromSmarts('O=C-O')

    for ring in atom_rings:
        if 10 <= len(ring) <= 18:  # Expanded range to include broader macrolide sizes
            # Create a sub-molecule of just this ring to match patterns within it
            submol = Chem.PathToSubmol(mol, ring)
            if submol.HasSubstructMatch(lactone_pattern):
                return True, "Contains macrocyclic lactone ring with an ester linkage"
    
    return False, "No characteristic macrolide macrocyclic lactone ring found"