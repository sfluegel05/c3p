"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:16beta-hydroxy steroid
A 16-hydroxy steroid in which the hydroxy group at position 16 has a beta-configuration.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for stereochemistry
    mol = Chem.AddHs(mol)
    
    try:
        # More flexible steroid core patterns that allow for variations
        steroid_patterns = [
            # Basic 4-ring system with flexible connectivity
            "*~1~*~*~*~2~*~*~*~3~*~*~*~4~*~*~*~*~4~*~3~*~2~1",
            
            # Alternative pattern allowing for different bond types
            "*~1~*~*~*~2~*~*~*~3~*~*~*~4~*~*~*~*~4~*~3~*~2~1",
            
            # Pattern for modified steroids
            "*~1~*~*~*~2~*~*~*~3~*~*~*~*~*~3~*~2~1"
        ]
        
        has_steroid_core = False
        for pattern in steroid_patterns:
            steroid_core = Chem.MolFromSmarts(pattern)
            if steroid_core is not None and mol.HasSubstructMatch(steroid_core):
                has_steroid_core = True
                break
                
        if not has_steroid_core:
            return False, "No steroid-like ring system found"

        # Patterns for 16-beta-hydroxy group with different possible environments
        beta_oh_patterns = [
            # General pattern for 16β-OH with beta stereochemistry
            "[C,c]~1~[C,c]~[C,c]~[C,c]~2~[C,c]~[C,c]~[C,c]~3~[C@@H](O)~[C,c]~[C,c]~[C,c]~3~[C,c]~2~1",
            
            # Alternative pattern focusing on the D-ring with 16β-OH
            "[C,c]~1~[C,c]~[C,c]~2~[C@@H](O)~[C,c]~[C,c]~[C,c]~2~1",
            
            # More specific pattern for common steroid variations
            "[C,c]~1~[C,c]~[C,c]~[C,c](@[C,c])~2~[C@@H](O)~[C,c]~[C,c]~2~1"
        ]
        
        has_16beta_oh = False
        for pattern in beta_oh_patterns:
            beta_oh = Chem.MolFromSmarts(pattern)
            if beta_oh is not None and mol.HasSubstructMatch(beta_oh, useChirality=True):
                has_16beta_oh = True
                break
                
        if not has_16beta_oh:
            return False, "No 16-beta-hydroxy group found or incorrect stereochemistry"

        # Count carbons and check for reasonable steroid size
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count < 10:  # More permissive minimum carbon count
            return False, f"Too few carbons ({carbon_count}) for a steroid structure"

        # Verify presence of hydroxy group
        hydroxy_pattern = Chem.MolFromSmarts("[OH]")
        if not mol.HasSubstructMatch(hydroxy_pattern):
            return False, "No hydroxy group found"

        # Check for reasonable molecular size
        if mol.GetNumAtoms() < 15:  # Including hydrogens
            return False, "Molecule too small to be a steroid"

        return True, "Contains steroid-like core with 16-beta-hydroxy group"

    except Exception as e:
        return False, f"Error in structure analysis: {str(e)}"