"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies polychlorinated dibenzodioxins and related compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxins
    and related compounds based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_member, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Count halogens (Cl, Br)
    num_cl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Cl]')))
    num_br = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Br]')))
    total_halogens = num_cl + num_br
    
    if total_halogens < 1:
        return False, "No halogens found"

    # PCDD pattern - dibenzodioxin core
    pcdd_pattern = Chem.MolFromSmarts('c1ccc2Oc3ccccc3Oc2c1')
    
    # PCDF pattern - dibenzofuran core
    pcdf_pattern = Chem.MolFromSmarts('c1ccc2oc3ccccc3c2c1')
    
    # PCB/PBB pattern - biphenyl core
    pcb_pattern = Chem.MolFromSmarts('c1ccccc1-c1ccccc1')
    
    # Formicamycin-like pattern (complex tricyclic system with specific substitution)
    formicamycin_pattern = Chem.MolFromSmarts('C1CC2=CC=CC=C2C(=O)C1')
    
    # Check for core structures
    is_pcdd = mol.HasSubstructMatch(pcdd_pattern)
    is_pcdf = mol.HasSubstructMatch(pcdf_pattern)
    is_pcb = mol.HasSubstructMatch(pcb_pattern)
    is_formicamycin = mol.HasSubstructMatch(formicamycin_pattern)
    
    # Classify based on core structure and halogenation
    if is_pcdd and total_halogens >= 2:
        return True, "Polychlorinated dibenzodioxin structure"
        
    if is_pcdf and total_halogens >= 2:
        return True, "Polychlorinated dibenzofuran structure"
        
    if is_pcb and total_halogens >= 2:
        # Check if halogens are attached to the aromatic rings
        hal_arom_pattern = Chem.MolFromSmarts('[Cl,Br]-c1ccccc1')
        if mol.HasSubstructMatch(hal_arom_pattern):
            return True, "Polychlorinated/brominated biphenyl structure"
            
    if is_formicamycin and total_halogens >= 1:
        # Additional checks for formicamycin-like structures
        if mol.HasSubstructMatch(Chem.MolFromSmarts('[OH]')):
            return True, "Formicamycin-like structure with halogen substitution"
            
    # Check for ambigol-like structures (biaryl ethers with halogens)
    if mol.HasSubstructMatch(Chem.MolFromSmarts('c1ccccc1Oc1ccccc1')) and total_halogens >= 2:
        if mol.HasSubstructMatch(Chem.MolFromSmarts('[OH]')):
            return True, "Ambigol-like halogenated biaryl ether"
            
    # If none of the above patterns match
    return False, "Does not match any recognized core structure pattern"