"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is defined as a C28 steroid with a modified side chain forming a lactone ring, and its substituted derivatives.

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
    
    # Define a more accurate steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts('C1CC2CCC3C(C2)C1CCC4=C3C=CC=C4')
    if steroid_pattern is None:
        return False, "Steroid pattern definition error"
    
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Define an improved lactone ring pattern (five-membered lactone ring)
    lactone_pattern = Chem.MolFromSmarts('O=C1COC(=O)C1')
    if lactone_pattern is None:
        return False, "Lactone pattern definition error"
    
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"
    
    # Ensure approximately C28 framework by counting carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 28:
        return False, f"Too few carbons for a C28 steroid, found {carbon_count}"
    
    return True, "Molecule has a steroid backbone, a lactone ring, and a C28 framework"