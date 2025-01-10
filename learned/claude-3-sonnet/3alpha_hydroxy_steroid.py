"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: 3alpha-hydroxy steroid
A 3-hydroxy steroid in which the 3-hydroxy substituent is in the alpha-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES with stereochemistry
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern (four fused rings)
    # More permissive pattern that matches both 5-alpha and 5-beta steroids
    steroid_core = Chem.MolFromSmarts(
        "[cH0,CH2,CH]1~[cH0,CH2,CH]~[cH0,CH2,CH]~[cH0,CH2,CH]2~[cH0,CH2,CH]~[cH0,CH2,CH]~[cH0,CH2,CH]3~[cH0,CH2,CH]~[cH0,CH2,CH]~[cH0,CH2,CH]4~[cH0,CH2,CH]~[cH0,CH2,CH]~[cH0,CH2,CH]~[cH0,CH2,CH]4~[cH0,CH2,CH]~3~[cH0,CH2,CH]~2~1"
    )
    
    if steroid_core is None:
        return False, "Error in steroid core SMARTS pattern"
        
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"

    # Check for 3-alpha-hydroxy group
    # [H] for explicit hydrogen, @ for specific stereochemistry
    alpha_3_hydroxy = Chem.MolFromSmarts('[OH][C@H]1[CH2][CH2]')
    
    if alpha_3_hydroxy is None:
        return False, "Error in hydroxy group SMARTS pattern"
        
    # Use useChirality=True to enforce stereochemistry matching
    if not mol.HasSubstructMatch(alpha_3_hydroxy, useChirality=True):
        return False, "No 3-alpha hydroxy group found"

    # Additional structure validation
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Check for characteristic steroid ring sizes (should find 6,6,6,5 or similar)
    ring_sizes = sorted([len(r) for r in ring_info.AtomRings()])
    if not (5 in ring_sizes and ring_sizes.count(6) >= 2):
        return False, "Ring sizes not characteristic of steroid structure"

    # Count carbons (steroids typically have 17+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, "Too few carbons for steroid structure"

    # Verify the molecule has appropriate number of rings and is not too small or large
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if not (250 < mol_weight < 1000):
        return False, "Molecular weight outside typical steroid range"

    return True, "Contains steroid core with 3-alpha hydroxy group"