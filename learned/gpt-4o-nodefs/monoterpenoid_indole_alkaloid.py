"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a monoterpenoid indole alkaloid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enhanced SMARTS patterns for indole moiety to capture common structural motifs.
    indole_patterns = [
        Chem.MolFromSmarts("c1c[nH]c2cccc3c2c1ccc3"),  # Standard indole
        Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3c2cccc3"), # Substituted indole variants
        Chem.MolFromSmarts("c1cc[nH]c2c1cccc2"),        # Possible extended scaffolds
        Chem.MolFromSmarts("c1cnc2c(c1)cccc2")          # More variations of substitutions
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in indole_patterns):
        return False, "No indole moiety found, or misidentified"
    
    # Define a more general SMARTS pattern for the monoterpenoid-like moiety.
    isoprene_like_patterns = [
        Chem.MolFromSmarts("C=C(C)C"),          # Basic isoprene unit
        Chem.MolFromSmarts("C(C)(C)C=C"),       # Alternate representation
        Chem.MolFromSmarts("C=C(C)C(C)C")       # Common branched version
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in isoprene_like_patterns):
        return False, "No monoterpenoid features found"
    
    # Check for functional groups such as methoxy, ester, ether
    functional_groups_patterns = [
        Chem.MolFromSmarts("[OX2]C"),           # Methoxy
        Chem.MolFromSmarts("C(=O)OC"),          # Ester
        Chem.MolFromSmarts("C-O"),              # Ether groups
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in functional_groups_patterns):
        return False, "Lacks common functional groups of monoterpenoid indole alkaloids"

    return True, "Contains indole moiety and monoterpenoid features with relevant functional groups"