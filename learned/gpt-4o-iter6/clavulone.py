"""
Classifies: CHEBI:36092 clavulone
"""
from rdkit import Chem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    A clavulone is a class of esterified prostanoids obtained from marine corals.
    
    It's characterized by complex cyclic structures, multiple ester groups, conjugated diene
    arrangements, and may contain halogens. Commonly contains at least one cyclopentenone moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect sufficient ester groups
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Insufficient ester groups found"
    
    # Search for typical cyclic frameworks, such as a cyclopentenone group
    cyclopentenone_pattern = Chem.MolFromSmarts("C1([OH])=CCCC1=O")
    if not mol.HasSubstructMatch(cyclopentenone_pattern):
        return False, "No cyclopentenone or similar cycle found"

    # Look for potential conjugated systems that include dienes
    diene_pattern = Chem.MolFromSmarts("C=CC=CC=C")
    has_diene = mol.HasSubstructMatch(diene_pattern)
    if not has_diene:
        return False, "No suitable conjugated diene system found"

    # Check for halogen atoms as typical but non-mandatory descriptors
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I]")
    has_halogen = mol.HasSubstructMatch(halogen_pattern)

    # Classify with reasoning based on the presence of halogens
    if has_halogen:
        return True, "Matches clavulone structure with halogens"
    else:
        return True, "Matches clavulone structure potentially without halogens"