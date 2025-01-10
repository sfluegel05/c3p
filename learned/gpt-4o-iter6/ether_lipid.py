"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid has a glycerol-like backbone with one or more carbon atoms bonded to alkyl chains via ether linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced pattern for glycerol-like backbone allowing potential ether links
    # Targets a three-carbon backbone with attachment points for ether links
    backbone_pattern = Chem.MolFromSmarts("O[C@@H]([C@H](O)CO)CO")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No suitable glycerol-like backbone found"

    # General ether linkage pattern for detecting various ether connections
    # Also accounts for potential single bond variants with alkoxy groups
    ether_pattern = Chem.MolFromSmarts("[CX4,O]O[CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if not ether_matches:
        return False, "No ether linkages detected"

    # Check for phosphate groups (optional) 
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)

    if has_phosphate:
        return True, "Contains ether linkage with glycerol-like backbone and phosphate group"
    return True, "Contains ether linkage with glycerol-like backbone, but no phosphate group found"

# Testing the function with a sample SMILES string
# result, reason = is_ether_lipid("CCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(C)=O")
# print(result, reason)