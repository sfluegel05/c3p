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

    # Define core structure patterns with more flexible matching
    # Allow for substituted carbons in the aromatic rings
    dibenzodioxin_pattern = Chem.MolFromSmarts('c1c(*)c(*)c2Oc3c(*)c(*)c(*)c(*)c3Oc2c1')  
    dibenzofuran_pattern = Chem.MolFromSmarts('c1c(*)c(*)c2oc3c(*)c(*)c(*)c(*)c3c2c1')    
    biphenyl_pattern = Chem.MolFromSmarts('c1c(*)c(*)c(*)c(c1)-c1c(*)c(*)c(*)c(*)c1')     
    
    # Check if SMARTS patterns were created successfully
    if None in (dibenzodioxin_pattern, dibenzofuran_pattern, biphenyl_pattern):
        return False, "Error in SMARTS patterns"

    # Count halogens
    num_cl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Cl]')))
    num_br = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Br]')))
    total_halogens = num_cl + num_br
    
    if total_halogens == 0:
        return False, "No halogens present"

    # Check for core structures
    is_dibenzodioxin = mol.HasSubstructMatch(dibenzodioxin_pattern)
    is_dibenzofuran = mol.HasSubstructMatch(dibenzofuran_pattern)
    is_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)
    
    # Count rings and check aromaticity
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    if aromatic_rings < 2:
        return False, "Insufficient aromatic rings"

    # Count other substituents that might indicate natural products
    oh_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OH]')))
    ome_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OMe]')))
    other_subst = oh_count + ome_count

    # Natural product check - more permissive to allow some substitution
    # but still exclude highly functionalized molecules
    if other_subst > 4 and aromatic_rings > 3:
        return False, "Too many substituents and rings - likely a complex natural product"

    # Determine core type and check requirements
    if is_dibenzodioxin:
        if num_cl < 1 and num_br < 1:
            return False, "Insufficient halogenation for dibenzodioxin"
        core_type = "dibenzodioxin"
    elif is_dibenzofuran:
        if num_cl < 1 and num_br < 1:
            return False, "Insufficient halogenation for dibenzofuran"
        core_type = "dibenzofuran"
    elif is_biphenyl:
        if total_halogens < 2:
            return False, "Insufficient halogenation for biphenyl"
        core_type = "biphenyl"
    else:
        return False, "No characteristic core structure found"

    # Build classification reason
    halogen_desc = []
    if num_cl > 0:
        halogen_desc.append(f"{num_cl} chlorine")
    if num_br > 0:
        halogen_desc.append(f"{num_br} bromine")
    halogen_str = " and ".join(halogen_desc)
    
    reason = f"Contains {core_type} core with {halogen_str} substituents"
    
    return True, reason