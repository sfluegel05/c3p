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
    
    # Look for glycine-like moieties ([NX3]-[CH2]-[C(=O)]-[X])
    glycine_pattern = Chem.MolFromSmarts("[NX3][CH2]C(=O)[X]")
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    
    # Exclude molecules with more than one glycine-like moiety
    if len(glycine_matches) != 1:
        return False, "Incorrect number of glycine-like moieties"
    
    # Look for acyl group ([CX3](=O)[NX3])
    acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    
    # Check if the acyl group is connected to the glycine-like moiety
    for acyl_match in acyl_matches:
        acyl_carbon = acyl_match[0]
        glycine_nitrogen = glycine_matches[0][0]
        if mol.GetBondBetweenAtoms(acyl_carbon, glycine_nitrogen):
            # Additional checks for charged groups, unusual bonding, etc.
            # ...
            return True, "Contains an N-acylglycine moiety"
    
    return False, "Acyl group and glycine-like moiety not connected"