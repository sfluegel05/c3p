"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: nucleobase analogue
A molecule that can substitute for a normal nucleobase in nucleic acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic size checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, "Molecular weight too high for nucleobase analogue"
    
    num_heavy = mol.GetNumHeavyAtoms()
    if num_heavy > 40:  # Increased from 30 to allow for more modifications
        return False, "Too many heavy atoms for nucleobase analogue"
    
    # Count atoms - made ranges more permissive
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:  # Must have at least 2 nitrogens
        return False, "Too few nitrogen atoms for nucleobase analogue"
    
    # Ring analysis
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"
    
    # Core structure patterns - made more general
    nucleobase_cores = [
        # Pyrimidine cores (uracil, thymine, cytosine-like)
        # Allow for more variation in oxidation states and substitution
        "[n,N]1[c,C][c,C][n,N][c,C][c,C]1",  # Basic pyrimidine ring
        "[n,N]1[c,C][c,C][n,N][c,C]1",  # 5-membered variant
        
        # Purine cores (adenine, guanine-like)
        # More permissive patterns that allow for modifications
        "[n,N]1[c,C]2[n,N][c,C][n,N][c,C][c,C]2[n,N][c,C]1",  # Basic purine
        "[n,N]1[c,C]2[n,N][c,C][n,N][c,C]2[n,N][c,C]1",  # Modified purine
        "[n,N]1[c,C]2[n,N][n,N][c,C][c,C]2[n,N][c,C]1",  # Aza-variant
    ]
    
    # Check for core structures
    found_core = False
    matched_core = None
    for pattern in nucleobase_cores:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_core = True
            matched_core = pattern
            break
            
    if not found_core:
        return False, "No nucleobase core structure found"
    
    # Common modifications to look for
    modifications = [
        ("C(=O)", "carbonyl"),
        ("N", "amino"),
        ("O", "hydroxy"),
        ("S", "thio"),
        ("F,Cl,Br,I", "halogen"),
        ("C=C", "alkene"),
    ]
    
    found_mods = []
    for pattern, mod_name in modifications:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_mods.append(mod_name)
    
    # Must have at least one typical nucleobase modification
    if not found_mods:
        return False, "Missing characteristic nucleobase modifications"
    
    # Check for characteristic features
    characteristic_features = [
        ("[n,N]1[c,C]([c,C])[n,N][c,C](=O)", "pyrimidine with carbonyl"),
        ("[n,N]1[c,C]2[n,N][c,C][n,N][c,C](N)", "purine with amino"),
        ("[n,N]1[c,C]2[n,N][c,C][n,N][c,C](=O)", "purine with carbonyl"),
    ]
    
    found_features = []
    for pattern, feature_name in characteristic_features:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_features.append(feature_name)
    
    # Success case
    core_type = "purine-like" if len(matched_core) > 30 else "pyrimidine-like"
    mods_str = ", ".join(found_mods[:3])  # List first 3 modifications
    return True, f"Contains {core_type} core with {mods_str} modifications"