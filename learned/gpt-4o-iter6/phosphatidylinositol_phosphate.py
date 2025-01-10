"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate
    based on its SMILES string.
    
    A phosphatidylinositol phosphate is characterized by a glycerol backbone with fatty acids,
    an inositol ring (C6 sugar alcohol), and one or multiple phosphate groups attached to the inositol.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a glycerol backbone pattern with ester linkages to fatty acids
    # We need a flexible pattern that allows for variable ester chains
    glycerol_pattern = Chem.MolFromSmarts("C(COC(=O))OCC(=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with ester linkages found"

    # Check for inositol ring structure
    # This should match any form of an inositol ring, typically represented as a hexahydroxycyclohexane
    inositol_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check for at least one phosphate group attached
    phosphate_pattern = Chem.MolFromSmarts("PO(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphates attached to the molecule"

    # Verify significant phosphate connection to inositol ring
    # Accept a range of phosphate attachments from mono- to trisphosphate
    inositol_phosphate_pattern = Chem.MolFromSmarts("c1c(O)cc(OP(=O)(O)O)cc(OP(=O)(O)O)c1O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "No significant phosphate attachment on inositol"

    return True, "Molecule contains glycerol, phosphates, and an inositol core"