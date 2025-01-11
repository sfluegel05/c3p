"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
from rdkit import Chem

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine consists of a glycerol backbone, ester-linked fatty acids, 
    a phosphate group, and a quaternary amine (choline) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("O=P(O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for quaternary amine (choline)
    choline_pattern = Chem.MolFromSmarts("OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(choline_pattern):
        return False, "No choline group found"

    # Check for ester linkage indicative of fatty acids
    ester_linkage_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_linkage_pattern)
    if len(ester_matches) < 1:
        return False, "No ester-linked fatty acids found"

    return True, "Contains glycerol backbone, phosphate group, choline, and ester-linked fatty acids"