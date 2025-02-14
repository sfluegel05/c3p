"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: CHEBI:36134 N-acylglycine
An N-acyl-amino acid in which the amino acid specified is glycine.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycine moiety ([NH]-[CH2]-[COOH])
    glycine_pattern = Chem.MolFromSmarts("[NH][CH2][C](=O)[OX1,OX2]")
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    if not glycine_matches:
        return False, "Glycine moiety not found"
    
    # Look for acyl group ([R]-[C(=O)]-[N])
    acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3][CH2][C](=O)[OX1,OX2]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "Acyl group not found"
    
    # Check if the acyl group and glycine moiety are connected
    for acyl_match in acyl_matches:
        for glycine_match in glycine_matches:
            if acyl_match[1] == glycine_match[0]:
                return True, "Contains an N-acylglycine moiety"
    
    return False, "Acyl group and glycine moiety not connected"