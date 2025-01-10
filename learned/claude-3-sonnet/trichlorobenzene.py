"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Definition: Any member of the class of chlorobenzenes carrying three chloro substituents at unspecified positions.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule contains a trichlorobenzene substructure
    (benzene ring with exactly 3 chlorine substituents at any positions).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains a trichlorobenzene substructure, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"
    
    # Get all benzene rings
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    # Check each benzene ring for chlorine substituents
    for ring_atoms in benzene_matches:
        chlorine_count = 0
        ring_atom_set = set(ring_atoms)
        
        # Count chlorines attached to this ring
        for ring_atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(ring_atom_idx)
            # Check neighbors of each carbon in the ring
            for neighbor in atom.GetNeighbors():
                # Only count chlorines that aren't part of the ring atoms
                if neighbor.GetSymbol() == "Cl" and neighbor.GetIdx() not in ring_atom_set:
                    chlorine_count += 1
        
        # If we found exactly 3 chlorines on this ring
        if chlorine_count == 3:
            return True, "Contains benzene ring with exactly 3 chlorine substituents"
            
    return False, "No benzene ring with exactly 3 chlorine substituents found"