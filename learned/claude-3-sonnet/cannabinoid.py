"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count basic atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if num_carbons < 15:  # Increased minimum carbon requirement
        return False, "Too few carbons for cannabinoid structure"
    if num_oxygens < 1:
        return False, "Cannabinoids must contain oxygen"

    # Classical cannabinoid patterns
    thc_core = Chem.MolFromSmarts("c1c(O)cc(CCCCC)cc1") # Phenolic ring with pentyl chain
    cbd_core = Chem.MolFromSmarts("c1c(O)c(CC=C)cc(CCCCC)c1O") # CBD-like core
    cbg_core = Chem.MolFromSmarts("c1c(O)c(CC=CC)cc(CCCCC)c1O") # CBG-like core
    
    # Synthetic cannabinoid patterns
    jwh_core = Chem.MolFromSmarts("c1ccc2n(CCCCC)c(C(=O))cc2c1") # JWH-like core
    cp_core = Chem.MolFromSmarts("c1cc(O)c(C2CCC(O)CC2)cc1") # CP-like core
    
    # Endocannabinoid specific patterns
    anandamide_core = Chem.MolFromSmarts("CCCCC=CCC=CCC=CCC=CCCC(=O)NCCO") # Anandamide-like
    ag_core = Chem.MolFromSmarts("CCCCC=CCC=CCC=CCC=CCCC(=O)OC(CO)CO") # 2-AG-like

    # Check for characteristic substitution patterns
    has_thc_core = mol.HasSubstructMatch(thc_core)
    has_cbd_core = mol.HasSubstructMatch(cbd_core)
    has_cbg_core = mol.HasSubstructMatch(cbg_core)
    has_jwh_core = mol.HasSubstructMatch(jwh_core)
    has_cp_core = mol.HasSubstructMatch(cp_core)
    has_anandamide = mol.HasSubstructMatch(anandamide_core)
    has_ag = mol.HasSubstructMatch(ag_core)

    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    
    # Classical cannabinoids
    if has_thc_core or has_cbd_core or has_cbg_core:
        if 250 < mol_wt < 400 and ring_count >= 2:
            return True, "Contains classical cannabinoid core structure"
            
    # Synthetic cannabinoids
    if has_jwh_core or has_cp_core:
        if 300 < mol_wt < 450 and ring_count >= 2:
            return True, "Contains synthetic cannabinoid core structure"
            
    # Endocannabinoids
    if has_anandamide or has_ag:
        if 300 < mol_wt < 400:
            return True, "Contains endocannabinoid core structure"

    # Additional checks for complex cannabinoid derivatives
    complex_core = Chem.MolFromSmarts("c1c(O)c([CH2,CH])cc(CCCC)c1") # Simplified core
    if mol.HasSubstructMatch(complex_core):
        if 250 < mol_wt < 500 and ring_count >= 2:
            # Check for key cannabinoid features
            phenol_pattern = Chem.MolFromSmarts("[OH]c1ccccc1")
            alkyl_chain = Chem.MolFromSmarts("CCCC")
            if mol.HasSubstructMatch(phenol_pattern) and mol.HasSubstructMatch(alkyl_chain):
                return True, "Contains modified cannabinoid core with required structural features"

    return False, "Does not match cannabinoid structural patterns"