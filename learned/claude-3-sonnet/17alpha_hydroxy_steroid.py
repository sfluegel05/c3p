"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern - four fused rings
    steroid_core = Chem.MolFromSmarts(
        '[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~1'
    )
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Multiple SMARTS patterns for 17α-OH to catch different representations
    oh_17_patterns = [
        # Pattern 1: Basic 17-OH with alpha orientation
        '[C;R1]1([C;R1][C;R1][C;R1]2)([O;H1])[C;R1][C;R1][C;R1]2',
        # Pattern 2: Explicit H version
        '[C;R1]1([C;R1][C;R1][C;R1]2)([O;H1])[C;R1][C;R1][C;R1]2[H]',
        # Pattern 3: Alternative representation
        '[C;R1]([O;H1])([C;R1])([C;R1])[C;R1]'
    ]
    
    oh_17_found = False
    for pattern in oh_17_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            oh_17_found = True
            break
    
    if not oh_17_found:
        return False, "No 17-alpha hydroxyl group found"

    # Validation checks
    # Count carbons (steroids typically have 19+ carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 19:
        return False, "Too few carbons for steroid structure"

    # Count rings
    ri = mol.GetRingInfo()
    if ri.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Check molecular weight
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 1200:  # Increased upper limit for glycosides
        return False, "Molecular weight outside typical range for steroids"

    # Count oxygens (should have at least one for the 17-OH)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found"

    # Additional check for tetracyclic system characteristic of steroids
    ring_atoms = set()
    for ring in ri.AtomRings():
        ring_atoms.update(ring)
    if len(ring_atoms) < 16:  # Typical steroid has at least 16 atoms in ring system
        return False, "Ring system too small for steroid structure"

    # Check for reasonable number of sp3 carbons
    sp3_carbons = sum(1 for atom in mol.GetAtoms() 
                     if atom.GetAtomicNum() == 6 and 
                     atom.GetHybridization() == Chem.HybridizationType.SP3)
    if sp3_carbons < 6:
        return False, "Insufficient sp3 carbons for steroid structure"

    return True, "Contains steroid core with 17-alpha hydroxyl group"