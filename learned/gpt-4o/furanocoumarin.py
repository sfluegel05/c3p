"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: Furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coumarin basic pattern (lactone ring with benzene fused)
    coumarin_pattern = Chem.MolFromSmarts("O=C2c1ccccc1OC2")

    # General furan pattern
    furan_pattern = Chem.MolFromSmarts("c1occc1")

    # Extended pattern for a furanocoumarin (fused furan-benzofuranone)
    furanocoumarin_fusion_pattern = Chem.MolFromSmarts("c1oc2ccccc2-3c1oc(=O)cc3")

    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin structure found"
    
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"
    
    if not mol.HasSubstructMatch(furanocoumarin_fusion_pattern):
        return False, "No furanocoumarin fusion detected"

    return True, "Structure contains a furan ring fused with a coumarin"