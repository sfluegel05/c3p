"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a pyrrole unit
    pyrrole_pattern = Chem.MolFromSmarts("n1cccc1")
    
    # Find all matches of the pyrrole pattern in the molecule
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    
    # Extract unique ring systems containing pyrrole units
    ring_info = mol.GetRingInfo()
    pyrrole_ring_systems = []
    for ring in ring_info.AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(mol.GetAtomWithIdx(idx).GetIsAromaticHeteroatom() for idx in ring):
            pyrrole_ring_systems.append(set(ring))
    
    # Count the number of unique ring systems containing two or more pyrrole units
    polypyrrole_ring_systems = [ring_system for ring_system in pyrrole_ring_systems if sum(1 for atom_idx in ring_system if mol.GetAtomWithIdx(atom_idx).GetIsAromaticHeteroatom()) >= 2]
    
    # Classify as a polypyrrole if at least one qualifying ring system is found
    if polypyrrole_ring_systems:
        return True, f"Contains {len(polypyrrole_ring_systems)} ring system(s) with two or more pyrrole units"
    else:
        return False, "No ring systems with two or more pyrrole units found"