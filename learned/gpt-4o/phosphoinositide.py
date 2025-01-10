"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern with at least two long carbon chains (C18+)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)COP")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    acyl_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCC")
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_matches) < 2:
        return False, f"Fewer than 2 acyl chains found, got {len(acyl_matches)}"

    # Look for inositol ring pattern (6-membered ring with 5 hydroxy groups)
    inositol_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C(O)1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check for phosphorylation on the inositol hydroxy groups
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphorylation found on inositol"

    return True, "Contains glycerol backbone, inositol ring, and phosphorylations on inositol hydroxy groups"