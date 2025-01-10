"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene.
    Tetrachlorobenzene is a benzene ring with exactly 4 chlorine atoms attached.
    
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
        
    # Get basic properties
    num_atoms = mol.GetNumAtoms()
    if num_atoms > 10:  # C6 + 4Cl = 10 atoms for tetrachlorobenzene
        # Allow for a few additional atoms for substituted variants
        if num_atoms > 15:  # Definitely too large
            return False, "Molecule too large to be a tetrachlorobenzene"
    
    # Count total chlorines
    total_chlorines = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    if total_chlorines < 4:
        return False, "Too few chlorine atoms"
    
    # Find benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    if not benzene_matches:
        return False, "No benzene ring found"
        
    # For each benzene ring, check if it's a valid tetrachlorobenzene
    for benzene_match in benzene_matches:
        benzene_atoms = set(benzene_match)
        chlorine_count = 0
        other_substituents = 0
        
        # Check each carbon in the benzene ring
        for atom_idx in benzene_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            
            # Check all neighbors of this carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 17:  # Chlorine
                    chlorine_count += 1
                elif neighbor.GetIdx() not in benzene_atoms:  # Other substituent
                    other_substituents += 1
                    
        # Valid tetrachlorobenzene conditions:
        # 1. Exactly 4 chlorines
        # 2. Limited other substituents (allow up to 2 for substituted variants)
        # 3. Not part of a larger aromatic system
        if (chlorine_count == 4 and 
            other_substituents <= 2 and 
            len(benzene_matches) == 1):
            
            # Additional check for fused ring systems
            fused_rings_pattern = Chem.MolFromSmarts("c12ccccc1cccc2")
            if not mol.HasSubstructMatch(fused_rings_pattern):
                return True, "Found valid tetrachlorobenzene structure"
                
    return False, "No valid tetrachlorobenzene structure found"