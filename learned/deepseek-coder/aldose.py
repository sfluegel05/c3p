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

    # Check for aldehyde group (C=O) in open-chain form
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Ensure the aldehyde is connected to a carbon chain with multiple hydroxyl groups
        aldehyde_match = mol.GetSubstructMatch(aldehyde_pattern)
        aldehyde_carbon = aldehyde_match[0]
        hydroxyl_count = 0
        for neighbor in mol.GetAtomWithIdx(aldehyde_carbon).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                for neighbor_neighbor in neighbor.GetNeighbors():
                    if neighbor_neighbor.GetAtomicNum() == 8 and neighbor_neighbor.GetDegree() == 1:  # Hydroxyl
                        hydroxyl_count += 1
        if hydroxyl_count >= 2:
            return True, "Contains aldehyde group with multiple hydroxyl groups"

    # Check for hemiacetal form (cyclic structure with oxygen in the ring)
    hemiacetal_pattern = Chem.MolFromSmarts("[OX2;R][CX4;R][CX3;R](=O)")
    if mol.HasSubstructMatch(hemiacetal_pattern):
        # Ensure the hemiacetal is part of a carbon chain with multiple hydroxyl groups
        hemiacetal_match = mol.GetSubstructMatch(hemiacetal_pattern)
        hemiacetal_carbon = hemiacetal_match[1]
        hydroxyl_count = 0
        for neighbor in mol.GetAtomWithIdx(hemiacetal_carbon).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                for neighbor_neighbor in neighbor.GetNeighbors():
                    if neighbor_neighbor.GetAtomicNum() == 8 and neighbor_neighbor.GetDegree() == 1:  # Hydroxyl
                        hydroxyl_count += 1
        if hydroxyl_count >= 2:
            return True, "Contains hemiacetal group with multiple hydroxyl groups"

    return False, "No aldehyde or hemiacetal group found with multiple hydroxyl groups"