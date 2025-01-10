"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol contains one or more isoprene units and an alcohol (OH) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification or exclusion
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Refined isoprene unit pattern [CH2C(Me)=CHCH2]
    isoprene_unit = Chem.MolFromSmarts("C(=C)C(C)C")
    num_isoprene_units = len(mol.GetSubstructMatches(isoprene_unit))

    if num_isoprene_units == 0:
        return False, "No isoprene units found"

    # Check for presence of alcohol (OH) or diphosphate (OPO3) groups
    oh_group = Chem.MolFromSmarts("[OX2H1]")
    has_oh_group = mol.HasSubstructMatch(oh_group)
    
    # It's crucial to recognize phosphates, as many prenols are phosphorylated
    phosphate_group = Chem.MolFromSmarts("P([OX1-])(=O)(O)O")
    has_phosphate_group = mol.HasSubstructMatch(phosphate_group)

    if not has_oh_group and not has_phosphate_group:
        return False, "No alcohol (OH) or phosphate group found"
    
    # Return True if both conditions of prenols are satisfied
    return True, f"Contains {num_isoprene_units} isoprene unit(s) and a valid group (OH/phosphate)"