"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is a glycerophospholipid having the polar alcohol inositol esterified
    to the phosphate group at the sn-3 position of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define patterns for key structural components
    glycerol_backbone_pattern = Chem.MolFromSmarts("C(CO)CO[P,O](=O)(OCC1)C1")
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1O")
    # Phosphate group next to the glycerol
    phosphate_inositol_pattern = Chem.MolFromSmarts("O[C@H](COP(=O)(O)O[C@H]1C(O)C(O)C(O)C(O)C1O)CO")

    # Check for the presence of glycerol backbone with phosphate linking to inositol
    if not mol.HasSubstructMatch(phosphate_inositol_pattern):
        return False, "No correct phosphate linkage to inositol found attached to glycerol"

    # Check for fatty acid chains (represented simplistically by C(=O)C)
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)C")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:  # Expect at least one fatty acid chain
        return False, "Missing fatty acid chains"

    return True, "Contains glycerophosphoinositol structure (glycerol backbone, attached phosphate group to inositol, inositol, and fatty acids)"