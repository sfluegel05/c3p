"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is defined as a C28 steroid lactone with modified side chains forming a lactone ring.

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
    
    # Check for steroid core (4 fused rings pattern)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4(C3)CCCC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core structure found"
    
    # Check for presence of lactone ring
    lactone_pattern = Chem.MolFromSmarts("O=C1OC=CC1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"
    
    # Check for C28 structure
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 28:
        return False, "Molecule has fewer than 28 carbons"
    
    return True, "Contains steroid core with lactone ring and sufficient carbon atoms"

# Example usage
example_smiles = "CC1=C(C)C(=O)O[C@H](C1)[C@](C)(O)[C@@]12O[C@@H]1C[C@H]1[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3CC[C@]21C"
is_withanolide(example_smiles)