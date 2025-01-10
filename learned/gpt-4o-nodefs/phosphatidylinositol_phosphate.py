"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    Phosphatidylinositol phosphates (PIPs) are lipids containing a glycerol backbone, 
    two fatty acids, a phosphate group, and an inositol ring that is phosphorylated at 
    least on one position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved glycerol backbone pattern typically with phosphate: allowing flexibility with fatty chain attachments
    glycerol_backbone_with_ester = Chem.MolFromSmarts("C(OC(=O)C)COC(=O)C")  # Allows ester linkages
    if not mol.HasSubstructMatch(glycerol_backbone_with_ester):
        return False, "No glycerol backbone with ester linkages found"
        
    # Check for an inositol ring with one or more phosphate attachments
    inositol_with_phosphate = Chem.MolFromSmarts("C1(COP(O)(O)=O)C(O)C(O)C(O)C(O)C1O")
    phosphate_count = len(mol.GetSubstructMatches(inositol_with_phosphate))
    if phosphate_count < 1:
        return False, "Inositol ring with at least one phosphate not found"

    return True, "Structure consistent with phosphatidylinositol phosphate"