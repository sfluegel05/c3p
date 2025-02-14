"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxin or related compound.
    These compounds include:
    - Polychlorinated dibenzodioxins
    - Polychlorinated dibenzofurans
    - Polychlorinated biphenyls
    - Polybrominated biphenyls
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a polychlorinated dibenzodioxin or related compound
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count atoms and check complexity
    num_atoms = mol.GetNumAtoms()
    if num_atoms > 50:  # Exclude large complex molecules
        return False, "Molecule too complex - likely a natural product"
    
    # Count halogens
    num_cl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Cl]')))
    num_br = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Br]')))
    total_halogens = num_cl + num_br
    
    # Define core structure patterns
    dibenzodioxin_pattern = Chem.MolFromSmarts('c1ccc2Oc3ccccc3Oc2c1')  # Dibenzodioxin core
    dibenzofuran_pattern = Chem.MolFromSmarts('c1ccc2oc3ccccc3c2c1')    # Dibenzofuran core
    biphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c1ccccc1')          # Biphenyl core
    
    # Check for core structures
    is_dibenzodioxin = mol.HasSubstructMatch(dibenzodioxin_pattern)
    is_dibenzofuran = mol.HasSubstructMatch(dibenzofuran_pattern)
    is_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)
    
    if not (is_dibenzodioxin or is_dibenzofuran or is_biphenyl):
        return False, "No dibenzodioxin, dibenzofuran, or biphenyl core structure found"

    # Count other substituents
    oh_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OH]')))
    ome_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OMe]')))
    other_subst = oh_count + ome_count
    
    # Specific requirements for each core type
    if is_dibenzodioxin:
        if num_cl < 2:
            return False, "Insufficient chlorination for dibenzodioxin"
        if other_subst > 2:
            return False, "Too many non-halogen substituents for PCDDs"
        core_type = "dibenzodioxin"
        
    elif is_dibenzofuran:
        if num_cl < 2:
            return False, "Insufficient chlorination for dibenzofuran"
        if other_subst > 2:
            return False, "Too many non-halogen substituents for PCDFs"
        core_type = "dibenzofuran"
        
    elif is_biphenyl:
        if total_halogens < 2:
            return False, "Insufficient halogenation for biphenyl"
        if num_cl + num_br < 2:
            return False, "Insufficient chlorination/bromination"
        if other_subst > 1:
            return False, "Too many non-halogen substituents for PCBs/PBBs"
        core_type = "biphenyl"
    
    # Check for natural product-like features
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 3:  # Core structure should have 2-3 rings only
        return False, "Too many rings - likely a natural product"
        
    # Build classification reason
    halogen_desc = []
    if num_cl > 0:
        halogen_desc.append(f"{num_cl} chlorine")
    if num_br > 0:
        halogen_desc.append(f"{num_br} bromine")
    halogen_str = " and ".join(halogen_desc)
    
    reason = f"Contains {core_type} core with {halogen_str} substituents"
    
    return True, reason