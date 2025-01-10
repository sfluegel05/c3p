"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is an amino acid (glycine) with an acyl group attached to the nitrogen atom.

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

    # Glycine core pattern with potential variations
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O") 
    # Check for primary amine bound to two or potentially non-carbon atoms in glycines
    if not mol.HasSubstructMatch(glycine_pattern):
        return False, "No glycine core structure found"
    
    # N-acyl linkage
    acyl_glycine_pattern = Chem.MolFromSmarts("C(=O)NCC(=O)O") 
    # Confirm acyl chain structure bound directly to amide
    if not mol.HasSubstructMatch(acyl_glycine_pattern):
        return False, "No acyl group linked to glycine nitrogen in structure"

    return True, "Contains characteristic N-acylglycine structure"