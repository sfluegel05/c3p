"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone
Definition: A class of esterified prostanoids obtained from marine corals.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core patterns characteristic of clavulones
    # Pattern 1: Cyclopentenone core with halogen or acetoxy
    core_pattern1 = Chem.MolFromSmarts("[Cl,Br,I,O][C]1[C](=O)[C][C][C]1") # Halogenated core
    core_pattern2 = Chem.MolFromSmarts("O[C]1[C](=O)[C][C][C]1") # Hydroxylated core
    
    # Pattern 2: Characteristic acetoxy group
    acetoxy = Chem.MolFromSmarts("OC(=O)C")
    
    # Pattern 3: Long unsaturated carbon chain
    chain_pattern = Chem.MolFromSmarts("C/C=C/C[C,O]")
    
    # Pattern 4: Conjugated system with ester
    conj_ester = Chem.MolFromSmarts("C=CC=C[CH]CC(=O)OC")
    
    # Check core structure
    has_core = mol.HasSubstructMatch(core_pattern1) or mol.HasSubstructMatch(core_pattern2)
    if not has_core:
        return False, "Missing characteristic cyclopentenone core"
    
    # Must have at least one acetoxy group
    acetoxy_count = len(mol.GetSubstructMatches(acetoxy))
    if acetoxy_count < 1:
        return False, "Missing acetoxy group"
    
    # Must have characteristic chain
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing characteristic unsaturated carbon chain"
        
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 18:
        return False, "Too few carbons for clavulone"
    if o_count < 3:
        return False, "Too few oxygens for clavulone"
        
    # Check for presence of halogens or multiple acetoxy groups
    hal_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9,17,35,53])
    if hal_count == 0 and acetoxy_count < 2:
        return False, "Must have either halogen or multiple acetoxy groups"
        
    # Additional structural features typical of clavulones
    if not (mol.HasSubstructMatch(conj_ester) or acetoxy_count >= 2):
        return False, "Missing characteristic conjugated system or multiple acetoxy groups"

    return True, "Contains characteristic clavulone core with appropriate substitution pattern"