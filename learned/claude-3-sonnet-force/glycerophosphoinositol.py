"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophosphoinositol(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is a glycerophospholipid with the polar alcohol inositol esterified
    to the phosphate group at the sn-3 position of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a glycerophosphoinositol, False otherwise
        str: Reason for the classification
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Check for the presence of an inositol ring
        inositol_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1O")
        inositol_match = mol.GetSubstructMatches(inositol_pattern)
        if not inositol_match:
            return False, "No inositol ring found"

        # Check for the presence of a glycerol backbone
        glycerol_pattern = Chem.MolFromSmarts("C(CO)C(O)CO")
        glycerol_match = mol.GetSubstructMatches(glycerol_pattern)
        if not glycerol_match:
            return False, "No glycerol backbone found"

        # Check for the presence of a phosphate group
        phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
        phosphate_match = mol.GetSubstructMatches(phosphate_pattern)
        if not phosphate_match:
            return False, "No phosphate group found"

        # Check if the inositol ring is esterified to the phosphate group at the sn-3 position
        ester_pattern = Chem.MolFromSmarts("OC(=O)C(O)C(O)CO[P@@](=O)(O)OC1C(O)C(O)C(O)C(O)C1O")
        ester_match = mol.GetSubstructMatches(ester_pattern)
        if not ester_match:
            return False, "Inositol ring not esterified to phosphate group at sn-3 position"

        return True, "Molecule matches the criteria for a glycerophosphoinositol"

    except Exception as e:
        return False, f"Error: {str(e)}"