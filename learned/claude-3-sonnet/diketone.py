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
    A diketone contains exactly two ketone groups (C=O where C is bonded to carbons).

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

    # Ketone pattern: Carbon double bonded to oxygen, where the carbon:
    # - is not part of a carboxylic acid (OH-C=O)
    # - is not part of an ester (O-C=O)
    # - is not part of an amide (N-C=O)
    # - is bonded to at least one carbon
    ketone_pattern = Chem.MolFromSmarts("""
        [CX3](=[OX1])
        [#6]
        ![$([CX3](=[OX1])[OX2H1])]
        ![$([CX3](=[OX1])[OX2][#6])]
        ![$([CX3](=[OX1])[NX3])]
    """)
    
    # Find all ketone matches
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Count unique ketone carbons
    ketone_carbons = set(match[0] for match in ketone_matches)
    
    if len(ketone_carbons) != 2:
        return False, f"Found {len(ketone_carbons)} ketone groups, need exactly 2"

    # Special patterns to recognize
    alpha_diketone = Chem.MolFromSmarts("[#6]C(=O)C(=O)[#6]")
    cyclic_13_diketone = Chem.MolFromSmarts("[#6]1[#6]C(=O)[#6][#6]C1=O")
    cyclic_14_diketone = Chem.MolFromSmarts("[#6]1[#6]C(=O)[#6][#6]C(=O)[#6]1")
    
    # Check for special patterns
    if mol.HasSubstructMatch(alpha_diketone):
        return True, "Contains alpha-diketone pattern"
    if mol.HasSubstructMatch(cyclic_13_diketone):
        return True, "Contains 1,3-diketone pattern"
    if mol.HasSubstructMatch(cyclic_14_diketone):
        return True, "Contains 1,4-diketone pattern"

    # Validate each ketone carbon
    for carbon_idx in ketone_carbons:
        atom = mol.GetAtomWithIdx(carbon_idx)
        carbon_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
        if carbon_neighbors == 0:
            return False, "Ketone carbon must be bonded to at least one carbon atom"

    return True, "Contains exactly two ketone groups"