"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    The classification focuses on arsenic within organic structural contexts.
    
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
    
    # Define SMARTS patterns that commonly represent organometalloidal arsenic compounds
    organometalloidal_patterns = [
        "[As](O)(O)=O",    # Arsenic with double-bonded oxygen and hydroxyl groups (as in arsenates)
        "[As](C)(C)(O)",   # Arsenic bonded to carbons and oxygen
        "[As]([O-])([O-])", # Arsenic with multiple oxygen anions
        "[c][As](C)",      # Arsenic bonded to an aromatic carbon and alkyl carbon
        "[As]1C=CC=C[CH]1", # Common structural context (arsolane ring)
        "[As]=O",          # Includes oxidized arsenic context
        "c1c([As])ccc(c1)[As]",  # Arsenic within conjugated systems
    ]
    
    # Check for any structural pattern indicative of organometalloidal compounds
    for pattern in organometalloidal_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, f"Pattern '{pattern}' matched, classified as organometalloidal"
    
    return False, "Presence of arsenic but not in organometalloidal context"