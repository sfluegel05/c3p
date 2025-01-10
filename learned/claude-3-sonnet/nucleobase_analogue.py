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
    if mol_wt > 400:  # Lowered from 500
        return False, "Molecular weight too high for nucleobase analogue"
    
    num_heavy = mol.GetNumHeavyAtoms()
    if num_heavy > 25:  # Lowered from 30
        return False, "Too many heavy atoms for nucleobase analogue"
    
    # Must contain rings
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"
    
    # Count nitrogens and carbons
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if n_count < 2 or n_count > 7:
        return False, "Incorrect number of nitrogen atoms for nucleobase"
    if c_count < 4 or c_count > 12:
        return False, "Incorrect number of carbon atoms for nucleobase"
    
    # Core structure patterns
    nucleobase_patterns = [
        # Pyrimidine cores (uracil, thymine, cytosine)
        "[nX2,NX3]1[cX3,CX3][cX3,CX3][nX2,NX3][cX3,CX3][cX3,CX3]1", # Basic pyrimidine
        "[NX3]1[CX3]=,:[CX3][NX3][CX3]=,:[CX3]1", # Reduced pyrimidine
        
        # Purine cores (adenine, guanine)
        "[nX2,NX3]1[cX3,CX3]2[nX2,NX3][cX3,CX3][nX2,NX3][cX3,CX3][cX3,CX3]2[nX2,NX3][cX3,CX3]1",
        "[NX3]1[CX3]=,:[CX3]2[NX3][CX3]=,:[NX3][CX3]=,:[CX3]2[NX3][CX3]1",
        
        # Modified cores
        "[nX2,NX3]1[cX3,CX3][cX3,CX3][nX2,NX3][cX3,CX3,nX2,NX3][cX3,CX3,nX2,NX3]1", # Modified pyrimidine
        "[nX2,NX3]1[cX3,CX3]2[nX2,NX3][cX3,CX3][nX2,NX3][cX3,CX3,nX2][cX3,CX3]2[nX2,NX3][cX3,CX3]1" # Modified purine
    ]
    
    found_core = False
    for pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_core = True
            break
            
    if not found_core:
        return False, "No nucleobase core structure found"
    
    # Essential functional groups
    essential_groups = [
        ("[CX3](=O)[NX3]", "amide"), # Amide
        ("[NX3;H2,H1;!$(NC=O)]", "amine"), # Primary/secondary amine
        ("[#6]=O", "carbonyl"), # Carbonyl
        ("[NX2]=[CX3]", "imine"), # Imine
    ]
    
    found_groups = []
    for pattern, group_name in essential_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_groups.append(group_name)
    
    if not found_groups:
        return False, "Missing characteristic nucleobase functional groups"
    
    # Additional filters
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)
    
    # Most atoms should be in rings for nucleobases
    ring_atom_ratio = len(ring_atoms) / num_heavy
    if ring_atom_ratio < 0.5:
        return False, "Too few atoms in ring systems"
    
    # Success case
    reason = f"Contains nucleobase core structure with {', '.join(found_groups)}"
    return True, reason