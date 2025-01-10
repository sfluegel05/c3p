"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for O-acyl L-carnitine core structure with key patterns
    # 1. Quaternary ammonium group [N+](C)(C)C
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+](C)(C)C")
    if not mol.HasSubstructMatch(quaternary_ammonium_pattern):
        return False, "No quaternary ammonium ion pattern found"

    # 2. Ester linkage pattern typical of acyl-carnitine
    ester_linkage_pattern = Chem.MolFromSmarts("OC(=O)C")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage pattern found"

    # Check for L-configuration at specific chiral centers
    expected_l_chiral_centers = [(atom.GetIdx(), atom.GetChiralTag()) 
                                 for atom in mol.GetAtoms() 
                                 if atom.GetChiralTag() == Chem.CHI_TETRAHEDRAL_CCW]
    
    if not expected_l_chiral_centers:
        return False, "No chiral centers with L-configuration found"

    # Additional checks could include:
    # - Expected bond patterns and lengths consistent with O-acyl-L-carnitine molecules
    # - Allowed functional groups without additional substitutions typical of other classes

    return True, "Contains key features of O-acyl-L-carnitine, including quaternary ammonium ion, ester linkage, and L-configuration"