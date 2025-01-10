"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    We focus on identifying the presence of arsenic bonded within typical organic structures.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of Arsenic (As)
    has_arsenic = any(atom.GetAtomicNum() == 33 for atom in mol.GetAtoms())
    if not has_arsenic:
        return False, "No arsenic atom found, not an organometalloidal compound"
    
    # Define SMARTS patterns for typical organometalloidal features
    organometalloidal_patterns = [
        "[As]=O",    # Arsenic with double bonded oxygen, a common feature
        "[As][C]",   # Arsenic directly bonded to carbon
        "[As][C][C]", # Extended bonding with carbons to form a chain
        "[As][C][N]", # Arsenic with carbon and heteroatom (organometallic environment)
        "[As][C][O]", # Functional groups like arsenates
        "c[As]",     # Aromatic carbon (phenyl) attached to arsenic
    ]
    
    # Check for any structural pattern indicative of organometalloidal compounds
    for pattern in organometalloidal_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, f"Pattern '{pattern}' matched, classified as organometalloidal"
    
    return False, "Arsenic present but not in an organometalloidal context"

# Example Usage:
# print(is_organometalloidal_compound("C[As](O)(O)=O"))  # Should return True with a valid reason