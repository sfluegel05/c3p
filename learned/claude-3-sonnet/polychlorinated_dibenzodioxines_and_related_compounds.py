"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies polychlorinated dibenzodioxins and related compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, "Insufficient ring count"
        
    # Count atoms and get basic properties
    num_atoms = mol.GetNumAtoms()
    num_carbons = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    num_cl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Cl]')))
    num_br = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Br]')))
    total_halogens = num_cl + num_br
    
    if total_halogens < 2:
        return False, "Insufficient halogenation"

    # Calculate halogen to carbon ratio
    halogen_carbon_ratio = total_halogens / num_carbons
    if halogen_carbon_ratio > 1.0:
        return False, "Too many halogens relative to carbons"

    # More specific patterns
    pcdd_pattern = Chem.MolFromSmarts('c1ccc2Oc3c(Cl,Br)cc(Cl,Br)cc3Oc2c1')
    pcdf_pattern = Chem.MolFromSmarts('c1ccc2oc3c(Cl,Br)cc(Cl,Br)cc3c2c1')
    pcb_pattern = Chem.MolFromSmarts('c1ccc(-c2ccc(Cl,Br)cc2)cc1')
    
    # Check for specific exclusion patterns (common false positives)
    glycopeptide_pattern = Chem.MolFromSmarts('NC(=O)C(N)C')
    if mol.HasSubstructMatch(glycopeptide_pattern) and num_atoms > 50:
        return False, "Appears to be a glycopeptide"

    # Check core structures with more specific requirements
    if mol.HasSubstructMatch(pcdd_pattern):
        # Verify proper halogen positions for PCDD
        if len(mol.GetSubstructMatches(Chem.MolFromSmarts('c1c(Cl,Br)cc2Oc3c(Cl,Br)cccc3Oc2c1'))) > 0:
            return True, "Polychlorinated dibenzodioxin structure"
            
    if mol.HasSubstructMatch(pcdf_pattern):
        # Verify proper halogen positions for PCDF
        if len(mol.GetSubstructMatches(Chem.MolFromSmarts('c1c(Cl,Br)cc2oc3c(Cl,Br)cccc3c2c1'))) > 0:
            return True, "Polychlorinated dibenzofuran structure"
            
    # More specific PCB/PBB pattern check
    if mol.HasSubstructMatch(pcb_pattern):
        # Additional checks for proper PCB/PBB structure
        aromatic_rings = len(mol.GetSubstructMatches(Chem.MolFromSmarts('c1ccccc1')))
        if aromatic_rings == 2 and total_halogens >= 2:
            # Check if most substituents are halogens
            non_hal_subst = len(mol.GetSubstructMatches(Chem.MolFromSmarts('c1ccccc1[!C;!Cl;!Br;!H]')))
            if non_hal_subst <= 2:  # Allow limited non-halogen substituents
                return True, "Polychlorinated/brominated biphenyl structure"
    
    # Count oxygen atoms in aromatic rings
    arom_oxy = len(mol.GetSubstructMatches(Chem.MolFromSmarts('c1ccccc1O')))
    if arom_oxy > 2 and ring_count > 3:
        return False, "Complex oxygenated aromatic structure - not a target compound"

    return False, "Does not match target compound patterns"