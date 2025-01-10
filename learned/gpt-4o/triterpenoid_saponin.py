"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define expanded triterpenoid backbone pattern
    # Example considering a broader set of possible sterane frameworks
    triterpenoid_patterns = [
        Chem.MolFromSmarts("C1CCC2(C)C3CCC4C(C)C(C)CCC4(C)C3CCC2C1"),  # protostane pattern
        Chem.MolFromSmarts("C1=CCC2(C)C(C)C3CCC4C(C)CCC4CCC3(C)C2C1"),  # hopane-like
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in triterpenoid_patterns):
        return False, "No triterpenoid backbone found"
    
    # Define glycosidic linkage presence
    # Check if there is a ring structure with oxygen indicating a sugar component
    glycosidic_oxygen = Chem.MolFromSmarts("O-C-C-O")  # Simplified pattern for minimal sugar involvement
    
    if not mol.HasSubstructMatch(glycosidic_oxygen):
        return False, "No glycosidic linkage with sugar rings found"
    
    return True, "Contains triterpenoid structure with glycosidic linkage"