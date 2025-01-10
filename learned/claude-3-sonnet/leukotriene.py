"""
Classifies: CHEBI:25029 leukotriene
"""
"""
Classifies: CHEBI:6051 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for C20 backbone
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Not a C20 backbone structure"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bonds < 4:
        return False, f"Insufficient double bonds (found {double_bonds}, need at least 4)"
    
    # Look for conjugated triene system (three consecutive double bonds)
    conjugated_triene = Chem.MolFromSmarts("C=CC=CC=C")
    if not mol.HasSubstructMatch(conjugated_triene):
        return False, "No conjugated triene system found"
    
    # Additional checks for common leukotriene features
    
    # Check for hydroxyl groups (common in many leukotrienes)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    # Check for potential cysteinyl group (present in LTC4, LTD4, LTE4)
    cysteinyl_pattern = Chem.MolFromSmarts("SCC(N)C(=O)")
    has_cysteinyl = mol.HasSubstructMatch(cysteinyl_pattern)
    
    # Calculate the number of rotatable bonds to ensure flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Structure too rigid for leukotriene"
    
    # Verify molecular weight is in reasonable range for leukotrienes
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for leukotrienes"
    
    # Build classification reason
    features = []
    if hydroxyl_count > 0:
        features.append(f"{hydroxyl_count} hydroxyl groups")
    if has_cysteinyl:
        features.append("cysteinyl group")
    features_str = ", ".join(features)
    
    return True, f"C20 structure with conjugated triene system, carboxylic acid, and {features_str}"