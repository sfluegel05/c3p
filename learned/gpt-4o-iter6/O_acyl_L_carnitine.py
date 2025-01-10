"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    O-acyl-L-carnitine should have a quaternary ammonium, an ester linkage,
    a specific chiral L-configuration, and a carboxylate group adjacent to the chiral center.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for quaternary ammonium group [N+](C)(C)C
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+](C)(C)C")
    if not mol.HasSubstructMatch(quaternary_ammonium_pattern):
        return False, "No quaternary ammonium ion pattern found"
    
    # Check for ester linkage
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage pattern found"

    # Check for chiral center with L-stereochemistry
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    is_correct_stereo = any(stereo == 'S' for _, stereo in chiral_centers)
    
    if not is_correct_stereo:
        return False, "No correct L-stereochemistry found"

    # Check for carboxylate group adjacency to chiral center
    carboxyl_pattern = Chem.MolFromSmarts("C[C@@H](C(=O)[O-])C")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Carboxylate group not adjacent to chiral center"

    return True, "Contains O-acyl-L-carnitine features: quaternary ammonium, ester linkage, L-stereo, carboxylate."