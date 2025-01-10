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
    
    # Basic size/composition checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, "Molecular weight too high for nucleobase analogue"
    
    num_heavy = mol.GetNumHeavyAtoms()
    if num_heavy > 30:
        return False, "Too many heavy atoms for nucleobase analogue"
    
    # Count atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 2 or n_count > 7:
        return False, "Incorrect number of nitrogen atoms for nucleobase"
    if c_count < 2 or c_count > 12:
        return False, "Incorrect number of carbon atoms for nucleobase"
    if o_count > 4:
        return False, "Too many oxygen atoms for nucleobase"

    # Core structure patterns - more specific than before
    nucleobase_cores = [
        # Pyrimidine cores (uracil, thymine, cytosine-like)
        "[nX2,NX3]1[cX3,CX3][cX3,CX3][nX2,NX3]([cX3,CX3](=O))[cX3,CX3]1(=O)", # Uracil-like
        "[nX2,NX3]1[cX3,CX3][cX3,CX3][nX2,NX3][cX3,CX3](N)[cX3,CX3]1=O", # Cytosine-like
        
        # Purine cores (adenine, guanine-like)
        "[nX2,NX3]1[cX3,CX3]2[nX2,NX3][cX3,CX3][nX2,NX3][cX3,CX3](N)[cX3,CX3]2[nX2,NX3][cX3,CX3]1", # Adenine-like
        "[nX2,NX3]1[cX3,CX3]2[nX2,NX3][cX3,CX3][nX2,NX3][cX3,CX3](=O)[cX3,CX3]2[nX2,NX3][cX3,CX3]1", # Guanine-like
        
        # Modified cores
        "[nX2,NX3]1[cX3,CX3][cX3,CX3][nX2,NX3][cX3,CX3](=O)[nX2,NX3]1", # Modified pyrimidine
        "[nX2,NX3]1[nX2,NX3][cX3,CX3][nX2,NX3][cX3,CX3](=O)[cX3,CX3]1" # Azapyrimidine
    ]
    
    # Check for core structures
    found_core = False
    for pattern in nucleobase_cores:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_core = True
            break
            
    if not found_core:
        return False, "No nucleobase core structure found"
    
    # Required functional group patterns in specific positions
    essential_features = [
        # At least one of these must be present
        ("[nX2,NX3]1[cX3,CX3]([cX3,CX3])[nX2,NX3][cX3,CX3](=O)[cX3,CX3]1=O", "pyrimidine with carbonyls"),
        ("[nX2,NX3]1[cX3,CX3]2[nX2,NX3][cX3,CX3][nX2,NX3][cX3,CX3](N)[cX3,CX3]2[nX2,NX3][cX3,CX3]1", "purine with amine"),
        ("[nX2,NX3]1[cX3,CX3]2[nX2,NX3][cX3,CX3][nX2,NX3][cX3,CX3](=O)[cX3,CX3]2[nX2,NX3][cX3,CX3]1", "purine with carbonyl")
    ]
    
    found_features = []
    for pattern, feature_name in essential_features:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_features.append(feature_name)
    
    if not found_features:
        return False, "Missing characteristic nucleobase features"
    
    # Ring analysis
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"
        
    # Success case
    return True, f"Contains {found_features[0]} structure with appropriate modifications"