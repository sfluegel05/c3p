"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    It checks for features typical of saccharolipids: long fatty acid chains, diverse sugar moieties,
    and various ester/ether linkages.

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
    
    # Checking for long hydrocarbon chains (simplified long alkyl chain pattern)
    num_long_chains = sum(1 for atom in mol.GetAtoms()
                          if atom.GetAtomicNum() == 6 and atom.GetDegree() > 2)
    if num_long_chains < 12:  # Arbitrary threshold for long chains
        return False, "No long alkyl chains detected"
    
    # Broadly capture sugar moieties
    ring_info = mol.GetRingInfo()
    num_sugar_like_rings = sum(1 for ring in ring_info.BondRings()
                               if len(ring) == 5 or len(ring) == 6)  # 5 or 6-membered sugar rings
    if num_sugar_like_rings == 0:
        return False, "No sugar-like rings detected"
    
    # Check for ester and glycosidic linkages
    ester_glyco_linkage_patterns = [
        Chem.MolFromSmarts("[C](=O)[O][C]"),  # Ester
        Chem.MolFromSmarts("O[C@H]"),         # Glycosidic
        Chem.MolFromSmarts("[C](=O)O[C@@H]")  # Another common motif
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in ester_glyco_linkage_patterns):
        return False, "No ester or glycosidic linkages found"

    return True, "Contains typical saccharolipid features: long fatty acid chains, sugar moieties, and characteristic linkages"