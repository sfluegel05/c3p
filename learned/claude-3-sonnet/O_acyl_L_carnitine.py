"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine has an L-carnitine backbone with an acyl group attached via an ester bond.
    
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

    # Check for carnitine backbone with specific stereochemistry
    # [C@H] for L-carnitine (the central carbon)
    # Allow for deuterated methyl groups in the trimethylammonium
    l_carnitine_pattern = Chem.MolFromSmarts('[C@H]([CH2][C]([O-])=O)([CH2][N+]([CH3,CD3])([CH3,CD3])([CH3,CD3]))[OX2]')
    d_carnitine_pattern = Chem.MolFromSmarts('[C@@H]([CH2][C]([O-])=O)([CH2][N+]([CH3,CD3])([CH3,CD3])([CH3,CD3]))[OX2]')
    
    if mol.HasSubstructMatch(d_carnitine_pattern):
        return False, "Incorrect stereochemistry - D-carnitine found"
    
    if not mol.HasSubstructMatch(l_carnitine_pattern):
        return False, "Missing L-carnitine backbone structure"

    # Check for ester linkage
    # More specific pattern that requires the acyl group
    ester_pattern = Chem.MolFromSmarts('[C@H]([CH2][C]([O-])=O)([CH2][N+])O[C](=[O])[#6]')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Missing or incorrect ester linkage"

    # Verify basic structure
    # Count key functional groups
    carboxylate_pattern = Chem.MolFromSmarts('[C](=[O])[O-]')
    ester_pattern = Chem.MolFromSmarts('[C](=[O])[O][C]')
    quat_n_pattern = Chem.MolFromSmarts('[N+]([CH3,CD3])([CH3,CD3])([CH3,CD3])')
    
    n_carboxylate = len(mol.GetSubstructMatches(carboxylate_pattern))
    n_ester = len(mol.GetSubstructMatches(ester_pattern))
    n_quat_n = len(mol.GetSubstructMatches(quat_n_pattern))
    
    if n_carboxylate < 1:
        return False, "Missing carboxylate group"
    if n_ester < 1:
        return False, "Missing ester linkage"
    if n_quat_n != 1:
        return False, "Must have exactly one quaternary ammonium group"

    # Check overall charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 0:
        return False, f"Total charge must be 0, found {total_charge}"

    # Get all chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    # Verify stereochemistry of central carbon
    central_matches = mol.GetSubstructMatches(l_carnitine_pattern)
    if not central_matches:
        return False, "Cannot determine stereochemistry"
    
    central_carbon_idx = central_matches[0][0]
    for idx, config in chiral_centers:
        if idx == central_carbon_idx and config != 'R':  # R configuration gives L-carnitine
            return False, "Incorrect stereochemistry - must be L-configuration"

    return True, "Valid O-acyl-L-carnitine structure with correct L-configuration"