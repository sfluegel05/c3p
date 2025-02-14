"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxin or related compound.
    
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

    # Count halogens (Cl, Br)
    num_cl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Cl]')))
    num_br = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Br]')))
    total_halogens = num_cl + num_br
    
    if total_halogens < 2:
        return False, "Insufficient halogen atoms"

    # Check for aromatic rings
    if not mol.GetSubstructMatches(Chem.MolFromSmarts('c1ccccc1')):
        return False, "No aromatic rings found"

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
    
    # Determine the specific type
    core_type = []
    if is_dibenzodioxin:
        core_type.append("dibenzodioxin")
    if is_dibenzofuran:
        core_type.append("dibenzofuran")
    if is_biphenyl:
        core_type.append("biphenyl")
    
    # Build classification reason
    halogen_desc = []
    if num_cl > 0:
        halogen_desc.append(f"{num_cl} chlorine")
    if num_br > 0:
        halogen_desc.append(f"{num_br} bromine")
    halogen_str = " and ".join(halogen_desc)
    
    core_str = " or ".join(core_type)
    reason = f"Contains {core_str} core with {halogen_str} substituents"
    
    # Additional checks for specific subclasses
    formicamycin_pattern = Chem.MolFromSmarts('CC1(C)C2=C(O)C=C(O)C3=C2C(=O)[C@@]4(O)C(=O)C5=CC=C(OC)C(=C5C)[C@H]4C13')
    ambigol_pattern = Chem.MolFromSmarts('c1c(O)c(Cl)cc(Oc2ccccc2)c1')
    
    if mol.HasSubstructMatch(formicamycin_pattern):
        reason += " (Formicamycin class)"
    elif mol.HasSubstructMatch(ambigol_pattern):
        reason += " (Ambigol class)"

    return True, reason