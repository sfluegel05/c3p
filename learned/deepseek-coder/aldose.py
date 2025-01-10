"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is a polyhydroxy aldehyde (H[CH(OH)]nC(=O)H, n >= 2) or its intramolecular hemiacetal.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aldehyde group (C=O) in any form
    aldehyde_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Count total hydroxyl groups in molecule
        hydroxyl_count = sum(1 for atom in mol.GetAtoms() 
                           if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
        if hydroxyl_count >= 2:
            return True, "Contains aldehyde group with multiple hydroxyl groups"

    # Check for hemiacetal form (cyclic structure with oxygen in the ring)
    # More general pattern that matches any cyclic structure with an oxygen
    ring_pattern = Chem.MolFromSmarts("[OX2;R1]")
    if mol.HasSubstructMatch(ring_pattern):
        # Count total hydroxyl groups in molecule
        hydroxyl_count = sum(1 for atom in mol.GetAtoms() 
                           if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
        if hydroxyl_count >= 2:
            # Check if it's a sugar-like structure
            # Look for multiple -CH(OH)- groups
            sugar_pattern = Chem.MolFromSmarts("[CHX4]([OX2H1])")
            sugar_matches = mol.GetSubstructMatches(sugar_pattern)
            if len(sugar_matches) >= 2:
                return True, "Contains cyclic structure with multiple hydroxyl groups (hemiacetal form)"

    return False, "No aldehyde or hemiacetal group found with multiple hydroxyl groups"