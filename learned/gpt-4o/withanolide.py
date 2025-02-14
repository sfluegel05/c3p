"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    Withanolides are characterized by a steroidal skeleton with a lactone-containing side chain and C28 framework including optional modifications.

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
    
    # Steroid backbone, focusing on key rings with potential modifications
    steroid_pattern = Chem.MolFromSmarts('[#6]1[#6,#8][#6]2[#6]3[#6,#7][#6]4[#6](#[#6])C([#6])([#8,#7])[#6][#6,#7][#6](=[#6])O1')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Capture expanded lactone forming groups, considering cyclic structures
    lactone_pattern = Chem.MolFromSmarts('O=C1O[CX3][CX3](=[OX1])[#6][#6]1')
    if not mol.HasSubstructMatch(lactone_pattern) and not mol.HasSubstructMatch(Chem.MolFromSmarts('O=C1COC(=O)[CX3]1')):
        return False, "No lactone ring found"
    
    # Counting framework atoms expected in a special role-density formulation
    total_heavy_atoms = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1])
    if total_heavy_atoms < 28:
        return False, f"Too few heavy atoms for a C28 steroid framework, found {total_heavy_atoms}"
    
    return True, "Molecule has a steroid backbone, a lactone ring, and a C28-like framework"