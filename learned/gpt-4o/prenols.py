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
        return False, "Invalid SMILES string"

    # Check for the presence of an alcohol (OH) group
    oh_group = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(oh_group):
        return False, "No alcohol (OH) group found"

    # Check for isoprene units [CH2-C(Me)=CH-CH2]
    isoprene_unit = Chem.MolFromSmarts("C(=C[C])C")
    num_isoprene_units = len(mol.GetSubstructMatches(isoprene_unit))

    if num_isoprene_units == 0:
        return False, "No isoprene units found"
    
    # Return True if both alcohol group and isoprene units are found
    return True, f"Contains {num_isoprene_units} isoprene unit(s) and an alcohol (OH) group"