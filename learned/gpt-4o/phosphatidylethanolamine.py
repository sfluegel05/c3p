"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    Phosphatidylethanolamines have a glycerophospholipid structure with two hydrophobic chains,
    a glycerol backbone, a phosphate group, and an ethanolamine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for glycerol backbone pattern without chiral specifity
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)COP(O)(O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Check presence of phosphate group within the backbone: O=P(O)(O)
    phosphate_pattern = Chem.MolFromSmarts("COP(O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group connectivity not found"
    
    # Look for ethanolamine group, which is a key feature
    ethanolamine_pattern = Chem.MolFromSmarts("OCCN")
    if not mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "No ethanolamine group found"

    # Ensure presence of two ester linkages
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Expected 2 ester linkages, found {len(ester_matches)}"

    return True, "Molecule classified as phosphatidylethanolamine"