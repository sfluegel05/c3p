"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    It considers features typical of saccharolipids: long fatty acid chains, diverse sugar moieties,
    and characteristic linkages such as ester and glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Checking for long hydrocarbon chains
    chain_length_threshold = 16  # Minimum number of carbons
    chain_pattern = Chem.MolFromSmarts("C" + "~C" * (chain_length_threshold - 1))
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long alkyl chains detected"
    
    # Capture sugar moieties (common monosaccharides)
    sugar_patterns = [
        Chem.MolFromSmarts("C1([H][H])([O][C@H]([C@H]([C@H]([C@H]1[O]*)[O*])[O*])[O*])"),  # Pyranose structures
        Chem.MolFromSmarts("C1([H][H])([O][C@H]([C@H]([C@H]([C@H]1[O*])[O*])[O*])")   # Furanose structures
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns):
        return False, "No recognized sugar moieties detected"
    
    # Check for ester and glycosidic linkages
    linkage_patterns = [
        Chem.MolFromSmarts("[C](=O)[O][C]"),  # Ester
        Chem.MolFromSmarts("[O][C@H]([C@H][O])"),  # Glycosidic linkage
        Chem.MolFromSmarts("[C](=O)[O][C@@H]")  # Stereo-glycosidic motif
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in linkage_patterns):
        return False, "No ester or glycosidic linkages found"

    return True, "Contains typical saccharolipid features: long fatty acid chains, sugar moieties, and characteristic linkages"