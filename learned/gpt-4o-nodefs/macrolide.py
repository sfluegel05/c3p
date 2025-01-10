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

    # Get the molecular rings information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # SMARTS pattern for an ester within a ring. 
    # C=O where C and O are part of the same ring
    ester_pattern = Chem.MolFromSmarts('C(=O)O')

    for ring in atom_rings:
        if 12 <= len(ring) <= 16:  # Check typical macrolide ring sizes
            # Create a sub-molecule of just this ring to match patterns within it
            submol = Chem.PathToSubmol(mol, ring)
            if submol.HasSubstructMatch(ester_pattern):
                return True, "Contains macrocyclic lactone ring with an ester linkage"

    return False, "No characteristic macrolide macrocyclic lactone ring found"