"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: diketone
A compound containing exactly two ketone functionalities
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone contains exactly two ketone groups (C=O where C is bonded to two other carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific ketone pattern:
    # - Carbon double bonded to oxygen
    # - Carbon must have exactly 3 bonds (X3)
    # - Carbon must be connected to exactly two other carbons
    # - The carbons attached must not be part of C=O groups themselves
    ketone_pattern = Chem.MolFromSmarts("[#6;!$(C=O)][CX3](=[OX1])[#6;!$(C=O)]")
    
    # Patterns to exclude
    quinone_pattern = Chem.MolFromSmarts("C1(=O)C=CC(=O)C=C1")  # para-quinone
    orthoquinone_pattern = Chem.MolFromSmarts("C1(=O)C(=O)C=CC=C1")  # ortho-quinone
    alpha_diketone_pattern = Chem.MolFromSmarts("[#6]C(=O)C(=O)[#6]")  # alpha-diketone
    cyclic_diketone_pattern = Chem.MolFromSmarts("C1C(=O)CC(=O)C1")  # cyclic 1,3-diketone
    
    # Find all ketone matches
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Count unique ketone carbons (some patterns might match the same carbon twice)
    ketone_carbons = set(match[1] for match in ketone_matches)
    
    if len(ketone_carbons) != 2:
        return False, f"Found {len(ketone_carbons)} ketone groups, need exactly 2"

    # Check for problematic patterns
    if mol.HasSubstructMatch(quinone_pattern) or mol.HasSubstructMatch(orthoquinone_pattern):
        return False, "Contains quinone structure"

    # Get all atoms involved in ketone groups
    ketone_atoms = set()
    for match in ketone_matches:
        ketone_atoms.update(match)

    # Validate each ketone group
    for carbon_idx in ketone_carbons:
        atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Check that the carbon is not part of a ring
        if atom.IsInRing() and not mol.HasSubstructMatch(cyclic_diketone_pattern):
            neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            if neighbors.count(6) != 2:  # Must have exactly two carbons attached
                return False, "Ketone group is part of an invalid ring system"
        
        # Check the environment of each ketone carbon
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # If neighbor is carbon
                # Check that neighbor is not part of problematic functional groups
                if neighbor.GetIsAromatic():
                    return False, "Ketone is conjugated with aromatic system"
                if neighbor.IsInRing() and atom.IsInRing():
                    if not mol.HasSubstructMatch(cyclic_diketone_pattern):
                        return False, "Invalid cyclic ketone arrangement"

    # If we have alpha-diketone pattern, make sure it's valid
    if mol.HasSubstructMatch(alpha_diketone_pattern):
        matches = mol.GetSubstructMatches(alpha_diketone_pattern)
        for match in matches:
            if all(idx in ketone_carbons for idx in [match[1], match[2]]):
                return True, "Contains valid alpha-diketone pattern"

    # Additional check for 1,3-diketone pattern
    if mol.HasSubstructMatch(cyclic_diketone_pattern):
        return True, "Contains valid cyclic 1,3-diketone pattern"

    # If we've made it here, check that the ketones are well-separated
    return True, "Contains exactly two ketone groups"