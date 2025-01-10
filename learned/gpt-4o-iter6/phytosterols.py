"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants 
    and may vary in carbon side chains and/or presence or absence of a double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More inclusive steroid backbone pattern including common variations
    steroid_pattern = Chem.MolFromSmarts("C1C[C@H]2CC[C@@]3(C)C2CC[C@H]4[C@]3(CC=C[C@H]1)C(C)(C)O")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No inclusive tetracyclic steroid backbone found"
    
    # Check for specific hydroxyl group positions typical in phytosterols
    hydroxyl_pattern = Chem.MolFromSmarts("[C@@H]([C@H]3[C@]4([C@@H](C)CC[C@]4(C)CC3)O)[C@H](O)C")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups typically found on phytosterols"
    
    # Presence of multiple types of unsaturation preferred in side chains or rings
    ring_unsaturation_pattern = Chem.MolFromSmarts("C=C[C@]3(C)[C@@H](CCC=C4C4(C)C)CC3(C)C")
    side_chain_unsaturation_pattern = Chem.MolFromSmarts("\C=C\C")
    has_unsaturation = mol.HasSubstructMatch(ring_unsaturation_pattern) or mol.HasSubstructMatch(side_chain_unsaturation_pattern)
    
    # Conclude on elements of phytosterols
    if has_unsaturation:
        return True, "Contains elements of unsaturation and inclusive tetracyclic steroid backbone typical of phytosterols"
    
    return False, "Fails to meet typical phytosterol criteria despite structural similarity"