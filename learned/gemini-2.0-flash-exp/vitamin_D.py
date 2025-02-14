"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D based on its SMILES string.
    Vitamin D is a seco-steroid with a characteristic broken B-ring, triene system,
    and hydroxyl groups at specific positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 0. Check for steroid skeleton.
    steroid_skeleton = Chem.MolFromSmarts("[C;R5]1[C;R5][C;R6]2[C;R6]3[C;R6]([C;R6]2[C;R5]1)[C;R6][C;R5]3")
    if not mol.HasSubstructMatch(steroid_skeleton):
        return False, "Not a steroid structure (missing basic rings)"
    
    # 1. Check for the specific broken B-ring with triene system
    #   This pattern matches the characteristic ring system and the broken B ring.
    #   The specific triene system and the double bond within the B-ring break are important here.
    #   This pattern also includes the crucial carbons to match to a vitamin D.
    broken_b_ring_pattern = Chem.MolFromSmarts("[C;R6]1=[C;R6]-1[C;R5]~[C;R]~[C;R]=[C;R]~[C;R]=[C;R]~[C;R]")

    if not mol.HasSubstructMatch(broken_b_ring_pattern):
       return False, "Not a seco-steroid structure (missing broken B ring with triene)"

    # 2. Check for hydroxyl groups at typical positions
    # The alpha-hydroxyl group at position 1 is very important, and often 25 or 24
    hydroxyl_pattern_1 = Chem.MolFromSmarts("[C;R][OH1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern_1)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl group found (at least one is required)"

    hydroxyl_pattern_2 = Chem.MolFromSmarts("[C;R]~[C;R]~[C;R]~[C;R]([C;R])([OH1])") # hydroxyl at 25 or 24 carbon in side chain
    hydroxyl_matches_2 = mol.GetSubstructMatches(hydroxyl_pattern_2)

    if len(hydroxyl_matches) < 1 and len(hydroxyl_matches_2) < 1:
         return False, "No hydroxyl group found (at least one in ring or at side chain is required)"


    # 3. Check for characteristic side chain starting with -CH-CH2-CH2-CH-
    side_chain_pattern = Chem.MolFromSmarts("[CHX4,CHX3]~[CH2X4]~[CH2X4]~[CHX4,CHX3]")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if len(side_chain_matches) < 1:
         return False, "No characteristic side chain found"
    
    # 4. Molecular weight range check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 600:
        return False, f"Molecular weight out of range: {mol_wt}"
    
    # 5. Check number of Carbons and Oxygens (for D2/D3 variants)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if not (c_count in [27, 28]): #check D2 (28C) or D3 (27C)
        return False, f"Incorrect number of carbons: {c_count}"
    if o_count < 1:
        return False, "Must have at least one oxygen"

    return True, "Meets vitamin D structural criteria"