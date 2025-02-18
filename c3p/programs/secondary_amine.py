"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:32877 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the secondary amine pattern: [N]([C])([C])[H]
    # This pattern matches a nitrogen atom bonded to two carbon atoms and one hydrogen atom
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;H1]([CX4])([CX4])")
    
    # Check if the molecule contains the secondary amine pattern
    if mol.HasSubstructMatch(secondary_amine_pattern):
        # Exclude cases where the nitrogen is part of an amide, sulfonamide, imine, or other non-secondary amine functional groups
        amide_pattern = Chem.MolFromSmarts("[NX3]([CX3]=[OX1])")
        sulfonamide_pattern = Chem.MolFromSmarts("[NX3]([SX4](=[OX1])=[OX1])")
        imine_pattern = Chem.MolFromSmarts("[NX2]=[CX3]")
        enamine_pattern = Chem.MolFromSmarts("[NX3]([CX3]=[CX3])")
        
        if (mol.HasSubstructMatch(amide_pattern) or 
            mol.HasSubstructMatch(sulfonamide_pattern) or 
            mol.HasSubstructMatch(imine_pattern) or 
            mol.HasSubstructMatch(enamine_pattern)):
            return False, "Nitrogen is part of an amide, sulfonamide, imine, or enamine, not a secondary amine"
        
        return True, "Contains a nitrogen atom bonded to two carbon atoms and one hydrogen atom (secondary amine)"
    else:
        # Check for cases where the nitrogen is part of a ring or bonded to non-carbon atoms
        ring_nitrogen_pattern = Chem.MolFromSmarts("[NX3;H1]([CX4])([CX4])@*")
        if mol.HasSubstructMatch(ring_nitrogen_pattern):
            # Ensure the nitrogen in the ring is not part of an amide, sulfonamide, imine, or enamine
            if (mol.HasSubstructMatch(amide_pattern) or 
                mol.HasSubstructMatch(sulfonamide_pattern) or 
                mol.HasSubstructMatch(imine_pattern) or 
                mol.HasSubstructMatch(enamine_pattern)):
                return False, "Nitrogen in ring is part of an amide, sulfonamide, imine, or enamine, not a secondary amine"
            return True, "Contains a nitrogen atom bonded to two carbon atoms and one hydrogen atom in a ring (secondary amine)"
        
        return False, "No secondary amine pattern found"