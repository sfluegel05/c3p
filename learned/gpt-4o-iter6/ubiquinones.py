"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones are benzoquinones with methoxy groups and a polyisoprenoid side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a benzoquinone core with specific placements for oxygens
    benzoquinone_core_pattern = Chem.MolFromSmarts("COc1c(C)c(=O)c(C(=O)c1C=C)OC")
    if not mol.HasSubstructMatch(benzoquinone_core_pattern):
        return False, "No specifically substituted benzoquinone core moiety found"
        
    # Check for polyisoprenoid side chains (considering variability)
    # Attempting to capture single repeated units extending from core
    polyisoprenoid_patterns = [
        Chem.MolFromSmarts("C=C(C)CCC=C"),  # basic isoprene unit
        Chem.MolFromSmarts("C=C(C)C"),      # consider alterations/short patterns
    ]
    
    has_polyisoprenoid_side_chain = any(mol.HasSubstructMatch(pattern) for pattern in polyisoprenoid_patterns)
    if not has_polyisoprenoid_side_chain:
        return False, "No polyisoprenoid side chain found"
    
    return True, "Contains the specific characteristics of a ubiquinone: correctly structured benzoquinone core with methoxy groups, and a variable polyisoprenoid side chain"