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

    # Look for a three-carbon glycerol-like backbone with possible ether linkages
    backbone_pattern = Chem.MolFromSmarts("O[C@H](C)CO")
    
    # Check for modified backbone forms with varied substitution (e.g., any stereo centers etc.)
    alternative_backbone_pattern = Chem.MolFromSmarts("O[C]C(OC)CO")
    
    if not (mol.HasSubstructMatch(backbone_pattern) or mol.HasSubstructMatch(alternative_backbone_pattern)):
        return False, "No suitable glycerol backbone found (considering ether variations)"

    # Improved ether linkage pattern specific to R-O-R' with aliphatic chains
    ether_pattern = Chem.MolFromSmarts("[OX2]([CX4])[CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    if not ether_matches:
        return False, "No ether linkages detected within expected context"
        
    # Check if there is any phosphate group attached (commonly found, not mandatory)
    phosphate_pattern = Chem.MolFromSmarts("[OP](=[O])(O)O")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)

    if has_phosphate:
        return True, "Contains ether linkage with glycerol-like backbone and phosphate group"
    else:
        return True, "Contains ether linkage with glycerol-like backbone, but no phosphate group found"

# Testing the function with a sample SMILES string
# result, reason = is_ether_lipid("CCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(C)=O")
# print(result, reason)