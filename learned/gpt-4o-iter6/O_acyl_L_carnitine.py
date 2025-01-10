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

    # Check for quaternary ammonium group [N+](C)(C)C connected to an alpha carbon
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+]([C])[C]C")
    if not mol.HasSubstructMatch(quaternary_ammonium_pattern):
        return False, "No quaternary ammonium ion pattern found"
    
    # Check for ester linkage C(=O)O[C@H] indicating L-stereochemistry
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage pattern consistent with L-stereochemistry found"

    # Check for chiral centers specifically marked as L (R for D-configuration)
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    is_l_stereo = any(stereo in ['S', '?'] and atom_idx in [match[0] for match in mol.GetSubstructMatches(ester_pattern)] for atom_idx, stereo in chiral_centers)
    
    if not is_l_stereo:
        return False, "No correct L-stereochemistry found"

    # Check for a carboxylate group (C(=O)[O-]) following a charged nitrogen
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "Carboxylate group not found adjacent to chiral center"
    
    return True, "Contains O-acyl-L-carnitine features: quaternary ammonium, ester linkage with correct stereochemistry, and carboxylate."