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

    # Check for a generalized benzoquinone core
    benzoquinone_core_pattern = Chem.MolFromSmarts("O=C1C=CC(=O)C=C1")
    if not mol.HasSubstructMatch(benzoquinone_core_pattern):
        return False, "No benzoquinone core moiety found"
        
    # Check for presence of at least one methoxy group (O-C group)
    methoxy_pattern = Chem.MolFromSmarts("CO")
    if not mol.HasSubstructMatch(methoxy_pattern):
        return False, "No methoxy group found"

    # Check for presence of a polyisoprenoid side chain
    # Flexible enough to match a single isoprene unit minimum
    polyisoprenoid_pattern = Chem.MolFromSmarts("C=C(C)CCC=C")
    if not mol.HasSubstructMatch(polyisoprenoid_pattern):
        return False, "No polyisoprenoid side chain found"
    
    return True, "Matches the characteristics of a ubiquinone: contains a generalized benzoquinone core, methoxy group, and a polyisoprenoid side chain"