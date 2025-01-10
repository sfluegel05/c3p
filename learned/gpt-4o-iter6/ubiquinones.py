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

    # Check for a benzoquinone core with methoxy groups
    benzoquinone_core_with_methoxy_pattern = Chem.MolFromSmarts("O=C1C=C(C(=O)C=C1)OC")
    if not mol.HasSubstructMatch(benzoquinone_core_with_methoxy_pattern):
        return False, "No benzoquinone core with methoxy groups found"

    # Check for polyisoprenoid side chain
    # Here we allow up to 10 isoprenoid units to accommodate variability
    polyisoprenoid_pattern = Chem.MolFromSmarts("C=C(C)CCC=C(C)CCC(=C(C)CCC=C(C)CC)*")
    if not mol.HasSubstructMatch(polyisoprenoid_pattern):
        return False, "No polyisoprenoid side chain found"
    
    return True, "Contains the specific characteristics of a ubiquinone: benzoquinone core with methoxy groups and a polyisoprenoid side chain"