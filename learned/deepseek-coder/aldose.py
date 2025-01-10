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
        # Count hydroxyl groups (OH)
        hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
        if hydroxyl_count >= 2:
            return True, "Contains aldehyde group with multiple hydroxyl groups"
    
    # Check for hemiacetal form (cyclic structure with oxygen in the ring)
    hemiacetal_pattern = Chem.MolFromSmarts("[OX2;R][CX4;R][CX3;R](=O)")
    if mol.HasSubstructMatch(hemiacetal_pattern):
        # Count hydroxyl groups (OH)
        hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
        if hydroxyl_count >= 2:
            return True, "Contains hemiacetal group with multiple hydroxyl groups"
    
    # Broaden the hemiacetal pattern to include more variations
    broad_hemiacetal_pattern = Chem.MolFromSmarts("[OX2;R][CX4;R][CX3;R](=O)")
    if mol.HasSubstructMatch(broad_hemiacetal_pattern):
        # Count hydroxyl groups (OH)
        hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
        if hydroxyl_count >= 2:
            return True, "Contains broad hemiacetal group with multiple hydroxyl groups"
    
    # Check for cyclic structure with oxygen in the ring (hemiacetal form)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        # Ensure the ring contains an oxygen (hemiacetal form)
        ring_atoms = set()
        for ring in ring_info.AtomRings():
            ring_atoms.update(ring)
        oxygen_in_ring = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring_atoms)
        if oxygen_in_ring:
            # Count hydroxyl groups (OH)
            hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
            if hydroxyl_count >= 2:
                return True, "Contains cyclic structure with oxygen and multiple hydroxyl groups"

    return False, "No aldehyde or hemiacetal group found with multiple hydroxyl groups"