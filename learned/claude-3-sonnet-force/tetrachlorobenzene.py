"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: CHEBI:35741 tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is any member of the class of chlorobenzenes carrying
    four chloro groups at unspecified positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count chlorine atoms
    n_chlorine = sum(atom.GetAtomicNum() == 17 for atom in mol.GetAtoms())
    
    # Tetrachlorobenzene must have exactly 4 chlorine atoms
    if n_chlorine != 4:
        return False, f"Found {n_chlorine} chlorine atoms, need exactly 4"
    
    # Find ring systems
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Check if any ring system has 4 chlorine substituents
    for ring in rings:
        ring_atoms = set(ring)
        chlorine_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
        
        if len(ring_atoms.intersection(chlorine_atoms)) == 4:
            # Check if chlorine atoms are substituents on the ring
            substituents = [mol.GetAtomWithIdx(atom_idx).GetTotalNumHs(includeNeighbors=True) == 0 and mol.GetAtomWithIdx(atom_idx).GetImplicitValence() == 1 for atom_idx in ring_atoms.intersection(chlorine_atoms)]
            if all(substituents):
                return True, "Contains a ring system with 4 chlorine substituents"
    
    return False, "No ring system found with 4 chlorine substituents"