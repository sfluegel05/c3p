"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    Focuses on identifying arsenic in organic contexts.
    
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
        return False, "No arsenic atom found"
    
    # Define more specific SMARTS patterns for organometalloidal contexts
    organometalloidal_patterns = [
        "[As](C)(O)=O", # Organometalloidal arsenic with organic carbon linkage
        "[As](c)(O)=O", # Arsenic bonded to aromatic carbon and oxygen
        "[As](C)(C)(C)", # Tri-substituted arsenic with carbons
        "[As](C)(C)", # Arsenic primarily bonded to carbon chains
        "[c][As](C)", # Aromatic group bonded with arsenic and carbon
        "c1cc[cH]c[As]", # Benzene-like rings with arsenic substitution
    ]
    
    # Check for any structural pattern indicative of organometalloidal compounds
    for pattern in organometalloidal_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, f"Pattern '{pattern}' matched, classified as organometalloidal"
    
    return False, "Presence of arsenic but not in organometalloidal context"

# Example usage for testing:
# print(is_organometalloidal_compound("C[As](O)(O)=O"))  # Should return True