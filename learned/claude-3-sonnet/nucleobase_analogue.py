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
    
    # Check molecular weight - should be reasonable for a nucleobase analogue
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, "Molecular weight too high for nucleobase analogue"
    
    # Must contain at least one ring
    if not mol.GetRingInfo().NumRings():
        return False, "No rings found"
    
    # Count number of nitrogen atoms - nucleobases have multiple nitrogens
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, "Too few nitrogen atoms"
    
    # Look for basic nucleobase-like patterns
    patterns = [
        # Pyrimidine-like pattern (uracil, thymine, cytosine)
        "[#7]1[#6][#6][#7][#6][#6]1",
        # Purine-like pattern (adenine, guanine)
        "[#7]1[#6]2[#7][#6][#7][#6][#6]2[#7][#6]1",
        # Modified pyrimidine pattern
        "[#7]1[#6][#6][#7][#6,#7][#6]1",
        # Modified purine pattern
        "[#7]1[#6]2[#7][#6][#7][#6,#7][#6]2[#7][#6]1"
    ]
    
    found_base_pattern = False
    for pattern in patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_base_pattern = True
            break
            
    if not found_base_pattern:
        return False, "No nucleobase-like ring pattern found"
    
    # Look for common functional groups in nucleobases
    carbonyl_pattern = Chem.MolFromSmarts("[#6]=O")
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    imine_pattern = Chem.MolFromSmarts("[NX2]=[#6]")
    
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    has_amine = mol.HasSubstructMatch(amine_pattern)
    has_imine = mol.HasSubstructMatch(imine_pattern)
    
    if not (has_carbonyl or has_amine or has_imine):
        return False, "Missing characteristic nucleobase functional groups"
    
    # Count aromatic atoms
    num_aromatic = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if num_aromatic < 3:
        return False, "Insufficient aromatic character"
    
    # Check for reasonable number of substituents
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    if num_heavy_atoms > 30:
        return False, "Too many heavy atoms for nucleobase analogue"
        
    # Success case - provide detailed reason
    reason = "Contains nucleobase-like ring system with appropriate functional groups"
    if has_carbonyl:
        reason += ", carbonyl groups"
    if has_amine:
        reason += ", amine groups"
    if has_imine:
        reason += ", imine groups"
    
    return True, reason