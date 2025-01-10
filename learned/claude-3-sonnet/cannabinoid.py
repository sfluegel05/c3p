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

    # Classical cannabinoid patterns (more specific)
    thc_core = Chem.MolFromSmarts("c1c(O)cc(CCCCC)c2c1OC(C)(C)[C@H]1CC=C(C)[C@H]21") # THC core
    cbd_core = Chem.MolFromSmarts("c1c(O)c(CC2=CCCC(C(=C)C)C2)cc(CCCCC)c1O") # CBD core
    cbg_core = Chem.MolFromSmarts("c1c(O)c(CC=C(C)CCC=C(C)C)cc(CCCCC)c1O") # CBG core
    
    # Synthetic cannabinoid patterns
    indole_core = Chem.MolFromSmarts("c1ccc2c(c1)c(C(=O))cn2CCCCC") # Indole-based
    pyrrole_core = Chem.MolFromSmarts("c1ccc(CC(=O)c2[nH]c3ccccc3c2)cc1") # Pyrrole-based
    cp_core = Chem.MolFromSmarts("[OH]C1CCC([c,C]2ccc(O)cc2)CC1") # CP-like core
    
    # Endocannabinoid patterns
    fatty_amide = Chem.MolFromSmarts("CCCCC(=CC=CC=CC=CC=CC)CC(=O)NCCO") # Anandamide-like
    glycerol_ester = Chem.MolFromSmarts("OCC(O)COC(=O)CCCCC=CC=CC=CC=C") # 2-AG-like
    ether_core = Chem.MolFromSmarts("CCCCC=CC=CC=CC=CCCCOC(CO)CO") # 2-AG ether-like

    # Check for characteristic substitution patterns
    classical_cores = [
        (thc_core, "THC-like"), 
        (cbd_core, "CBD-like"),
        (cbg_core, "CBG-like")
    ]
    
    synthetic_cores = [
        (indole_core, "indole-based synthetic"),
        (pyrrole_core, "pyrrole-based synthetic"),
        (cp_core, "CP-like synthetic")
    ]
    
    endocannabinoid_cores = [
        (fatty_amide, "fatty acid ethanolamide"),
        (glycerol_ester, "monoacylglycerol"),
        (ether_core, "glyceryl ether")
    ]

    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Check for classical cannabinoids
    for core, core_type in classical_cores:
        if mol.HasSubstructMatch(core):
            if 280 < mol_wt < 400 and ring_count >= 2:
                return True, f"Classical cannabinoid with {core_type} core structure"

    # Check for synthetic cannabinoids
    for core, core_type in synthetic_cores:
        if mol.HasSubstructMatch(core):
            if 300 < mol_wt < 450 and ring_count >= 2:
                return True, f"Synthetic cannabinoid with {core_type} core structure"

    # Check for endocannabinoids
    for core, core_type in endocannabinoid_cores:
        if mol.HasSubstructMatch(core):
            if 300 < mol_wt < 500 and rotatable_bonds > 10:
                return True, f"Endocannabinoid ({core_type})"

    # Additional checks for MAGs and ethanolamides
    if "COC(=O)" in smiles and "CO" in smiles:  # Potential MAG
        if rotatable_bonds > 10 and "C=C" in smiles:  # Long unsaturated chain
            if 300 < mol_wt < 500:
                return True, "Monoacylglycerol (MAG) cannabinoid"
                
    if "NCCO" in smiles and "C(=O)" in smiles:  # Potential ethanolamide
        if rotatable_bonds > 10 and "C=C" in smiles:  # Long unsaturated chain
            if 300 < mol_wt < 500:
                return True, "Fatty acid ethanolamide cannabinoid"

    return False, "Does not match cannabinoid structural patterns"