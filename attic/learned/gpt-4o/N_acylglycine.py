"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine consists of a glycine moiety with an acyl group bonded to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycine core pattern (more flexible)
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)[O,N]")  # allows flexibility with acids or derivatives
    if not mol.HasSubstructMatch(glycine_pattern):
        return False, "No glycine core structure found"
    
    # N-acyl linkage pattern (specific to an acyl group attached via nitrogen)
    acyl_glycine_pattern = Chem.MolFromSmarts("C(=O)NCC(=O)[O,N]")  # carbonyl group bonded to nitrogen with glycine motif
    if not mol.HasSubstructMatch(acyl_glycine_pattern):
        return False, "No acyl group bonded to glycine nitrogen as expected for N-acylglycine"

    return True, "Contains N-acylglycine characteristic structure"