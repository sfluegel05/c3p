"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are terpenoids derived from sesquiterpenes with C15 skeletons,
    possibly rearranged or modified by the removal of one or more skeletal atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count total carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Fewer than 15 carbon atoms; not a typical sesquiterpenoid"
    
    # SMARTS pattern for core sesquiterpene structure (15 carbons in skeletal C-C bonds)
    core_pattern = Chem.MolFromSmarts("C1(C)CCCCC2(C=CCC=CC2)C1")
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing key sesquiterpene structural elements"
    
    # Checking for typical sesquiterpenoid functional groups 
    functional_patterns = [
        Chem.MolFromSmarts('O'),  # Hydroxyl groups
        Chem.MolFromSmarts('C=O'),  # Carbonyl groups (aldehydes/ketones)
        Chem.MolFromSmarts('O=C(O)'),  # Esters/lactones
    ]
    
    for pattern in functional_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a sesquiterpene skeleton with characteristic sesquiterpenoid functionalities"
    
    return False, "Does not possess common sesquiterpenoid functionalities"