"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: CHEBI:17596 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate consists of a nucleobase attached to a sugar with one or more phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of phosphorus (phosphate requirement)
    has_phosphorus = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
    if not has_phosphorus:
        return False, "No phosphorus atoms found"

    # Define nucleobase pattern (purine/pyrimidine core)
    base_pattern = Chem.MolFromSmarts("[n&R]1[c&R][n&R][c&R][n&R]1")  # Purine base
    if not mol.HasSubstructMatch(base_pattern):
        # Try pyrimidine pattern if purine not found
        base_pattern = Chem.MolFromSmarts("[n&R]1[c&R][n&R][c&R][c&R]1")  # Pyrimidine base
        if not mol.HasSubstructMatch(base_pattern):
            return False, "No nucleobase detected"

    # Define sugar pattern (ribose-like structure)
    sugar_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@@H](O)[C@H]1O")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar moiety detected"

    # Check phosphate linkage to sugar
    phosphate_attached = False
    for sugar_match in sugar_matches:
        # Get the oxygen atoms in the sugar
        sugar_o_indices = [i for i in sugar_match if mol.GetAtomWithIdx(i).GetAtomicNum() == 8]
        
        # Check if any sugar oxygen is connected to phosphorus
        for o_idx in sugar_o_indices:
            atom = mol.GetAtomWithIdx(o_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 15:  # Phosphorus
                    phosphate_attached = True
                    break
            if phosphate_attached:
                break
        if phosphate_attached:
            break

    if not phosphate_attached:
        return False, "No phosphate attached to sugar"

    return True, "Contains nucleobase-sugar with phosphate linkage"