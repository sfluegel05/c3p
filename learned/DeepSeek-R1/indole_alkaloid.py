"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: CHEBI:24828 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid must contain an indole skeleton (bicyclic structure with 
    benzene fused to pyrrole ring) and be an alkaloid (nitrogen-containing compound).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define indole core pattern (benzene fused to pyrrole)
    indole_pattern = Chem.MolFromSmarts("n1ccc2ccccc12")
    
    # Check for indole skeleton
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole skeleton found"

    # Basic alkaloid check: presence of at least one nitrogen
    # (indole's nitrogen already counts, but verify molecule is nitrogen-containing)
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 1:
        return False, "No nitrogen atoms found"

    # Additional check for basic amine nitrogen (not in amide or aromatic system)
    # This helps exclude simple indoles and non-alkaloid compounds
    amine_pattern = Chem.MolFromSmarts("[NX3;!$(N[a]);!$(N[C]=O)]")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No basic amine nitrogen detected"

    return True, "Contains indole skeleton with basic nitrogen characteristic of alkaloids"