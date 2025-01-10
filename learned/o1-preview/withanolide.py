"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is defined as any steroid lactone that is a C28 steroid 
    with a modified side chain forming a lactone ring and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define steroid nucleus SMARTS pattern (tetracyclic steroid backbone)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3(CCCC4)C')
    if steroid_pattern is None:
        return False, "Invalid steroid SMARTS pattern"
    
    # Check for steroid nucleus
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"
        
    # Define lactone ring SMARTS pattern (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts('C1OC(=O)[C;R]1')  # 5-membered lactone
    if lactone_pattern is None:
        return False, "Invalid lactone SMARTS pattern"
    
    # Find lactone rings
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone ring found"
    
    # Get the steroid nucleus atoms
    steroid_match = mol.GetSubstructMatch(steroid_pattern)
    steroid_atoms = set(steroid_match)
    
    # Check if any lactone ring is in the side chain (not part of steroid nucleus)
    side_chain_lactone = False
    for lactone_ring in lactone_matches:
        lactone_atoms = set(lactone_ring)
        if not lactone_atoms.issubset(steroid_atoms):
            side_chain_lactone = True
            break
    if not side_chain_lactone:
        return False, "Lactone ring is not in side chain"
    
    # Optional: check approximate carbon count (allowing some flexibility)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 26:
        return False, f"Too few carbons for a withanolide (found {c_count} carbons)"
    
    return True, "Molecule is a withanolide (steroid nucleus with side-chain lactone ring)"