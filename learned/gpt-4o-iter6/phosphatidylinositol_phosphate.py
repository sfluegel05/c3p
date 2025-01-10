"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate
    based on its SMILES string.
    
    A phosphatidylinositol phosphate is characterized by a glycerol backbone with esterified fatty acids,
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

    # Improved pattern for glycerol backbone with ester linkages
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with ester linkages found"

    # General pattern for an inositol ring (may include attached phosphates)
    inositol_phosphate_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1(OP(=O)(O)O)")
    matched = mol.GetSubstructMatches(inositol_phosphate_pattern)
    phosphate_attached_to_inositol = any(
        mol.mol.GetSubstructMatch(Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1OP(=O)(O)O"))
        for inositol_ring in matched)

    if not phosphate_attached_to_inositol:
        return False, "Inositol ring with phosphates not found"

    # Check for at least one ester linkage implying fatty acid
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if len(mol.GetSubstructMatches(ester_pattern)) < 2:
        return False, "Does not have at least two ester groups indicating fatty acids"

    return True, "Molecule is a phosphatidylinositol phosphate with appropriate glycerol, inositol, and phosphate structure"