"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    Phosphatidylinositols are characterized by a glycerol backbone linked to
    two fatty acids via ester bonds, a phosphate group, and a 1D-myo-inositol head.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a phosphatidylinositol, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the inositol ring structure - 1D-myo-inositol
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Look for the glycerol-like structure with linked phosphate group
    glycerol_phosphate_pattern = Chem.MolFromSmarts("C(COP(=O)(O)O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "Phosphate-linked glycerol backbone not found"

    # Check for at least two ester linkages (-C(=O)OC-)
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Less than two ester linkages found"

    return True, "Molecule matches phosphatidylinositol structure"

# Example execution
smiles_example = "P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)COC(=O)CCC/C=C\\C/C=C\\CCCC)(O)=O"
is_phosphatidylinositol(smiles_example)