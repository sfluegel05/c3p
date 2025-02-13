"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate
    based on its SMILES string.
    
    A phosphatidylinositol phosphate is characterized by a glycerol backbone with esterified fatty acids,
    an inositol ring (C6 sugar alcohol), and one or more phosphate groups attached to the inositol.
    
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

    # Pattern for glycerol backbone with ester linkages. Adjusting to ensure fatty acid tails.
    glycerol_esterified_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)C")
    if not mol.HasSubstructMatch(glycerol_esterified_pattern):
        return False, "No glycerol backbone with ester linkages found"

    # Pattern for an inositol ring with attached phosphate groups.
    inositol_phosphate_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1(OP(=O)(O)O)")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "Inositol ring with phosphates not found"

    # Check for at least two ester linkages indicating fatty acids
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if len(mol.GetSubstructMatches(ester_pattern)) < 2:
        return False, "Does not have at least two ester groups indicating fatty acids"

    return True, "Molecule is a phosphatidylinositol phosphate with appropriate glycerol, inositol, and phosphate structure"