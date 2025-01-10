"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    Biflavonoids generally contain two flavonoid moieties connected by covalent bonds.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General aromatic ring pattern and typical flavonoid backbone
    flavonoid_core_patterns = [
        Chem.MolFromSmarts("c1(ccc(c(c1)c2=cc(oc2=O)c)C(=O)C)"),  # example flavonoid core
        Chem.MolFromSmarts("c1cc(c2c(c1O)oc(=O)c(c2)C)"),  # another example pattern
        Chem.MolFromSmarts("C1=CC=C(C=C1)C2=CC(=O)C3=C(O2)C=CC=C3"), # basic flavan-3-ol
    ]

    # Check for at least two occurrences of flavonoid core patterns
    match_count = 0
    for core_pattern in flavonoid_core_patterns:
        matches = mol.GetSubstructMatches(core_pattern)
        match_count += len(matches)

    if match_count < 2:
        return False, "Less than two flavonoid units found"

    # Check for connections typical of biflavonoids, (C-O, C-C, or C=C between flavonoid units)
    linkage_patterns = [
        Chem.MolFromSmarts("[c]~[o]~[c]"),  # C-O-C linkage
        Chem.MolFromSmarts("[c]-[c]-[c]"),  # C-C linkage
        Chem.MolFromSmarts("[c]=[c]"),  # C=C linkage
    ]
    
    # Verify the presence of covalent linkage between flavonoid units
    has_linkage = any(mol.HasSubstructMatch(linkage) for linkage in linkage_patterns)

    if has_linkage:
        return True, "Contains dual flavonoid units connected by typical biflavonoid linkage patterns"
    
    return False, "Molecule does not consist of connected flavonoid units typical of biflavonoids"