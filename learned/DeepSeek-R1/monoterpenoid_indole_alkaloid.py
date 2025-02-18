"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: CHEBI:65312 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    The molecule should contain an indole moiety and structural features from a monoterpenoid (diisoprenoid).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for indole substructure (nH in the five-membered ring)
    indole_pattern = Chem.MolFromSmarts("[nH]1c2ccccc2cc1")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole moiety found"
    
    # Check for at least 20 carbons (approximate for monoterpenoid + indole)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Only {c_count} carbons, expected >=20"
    
    # Check for multiple rings (indole + terpenoid-derived structure)
    sssr = Chem.GetSSSR(mol)
    if len(sssr) < 2:
        return False, f"Only {len(sssr)} rings, expected multiple"
    
    # Check for at least one oxygen (common in terpenoid parts)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found"
    
    return True, "Contains indole moiety with terpenoid features"