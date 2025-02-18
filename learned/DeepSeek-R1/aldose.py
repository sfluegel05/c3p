"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: CHEBI:17476 aldose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is a polyhydroxy aldehyde (H[CH(OH)]nC(=O)H, n ≥ 2) or its cyclic hemiacetal form.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aldehydic form (open chain aldehyde with polyhydroxy chain)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1]=O")  # Aldehyde group
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Verify at least two hydroxyl groups in the molecule
        hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
        if hydroxyl_count >= 2:
            # Check for at least two CH(OH) groups adjacent to aldehyde
            adjacent_oh_count = 0
            for match in mol.GetSubstructMatches(aldehyde_pattern):
                ald_carbon = match[0]
                for neighbor in mol.GetAtomWithIdx(ald_carbon).GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        for ch_group in neighbor.GetNeighbors():
                            if ch_group.GetAtomicNum() == 8 and ch_group.GetTotalNumHs() >= 1:
                                adjacent_oh_count += 1
                                break
            if adjacent_oh_count >= 2:
                return True, "Aldehyde group with adjacent polyhydroxy chain found"
    
    # Check for cyclic hemiacetal form (oxygen in ring with adjacent hydroxyl)
    ring_hemiacetal_pattern = Chem.MolFromSmarts("[C;r][OH]-[O;r]")
    if mol.HasSubstructMatch(ring_hemiacetal_pattern):
        # Verify the hemiacetal carbon has only one non-ring neighbor (aldose characteristic)
        for match in mol.GetSubstructMatches(ring_hemiacetal_pattern):
            hemi_carbon = match[0]
            atom = mol.GetAtomWithIdx(hemi_carbon)
            non_ring_neighbors = [n for n in atom.GetNeighbors() if not n.IsInRing()]
            if len(non_ring_neighbors) == 1:  # Characteristic of aldose-derived hemiacetal
                return True, "Cyclic hemiacetal structure characteristic of aldose found"

    return False, "No aldehyde group or valid cyclic hemiacetal structure detected"