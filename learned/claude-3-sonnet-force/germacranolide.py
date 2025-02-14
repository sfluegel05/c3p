"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: CHEBI:36839 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on the germacrane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for germacrane skeleton patterns
    germacrane_patterns = [
        Chem.MolFromSmarts("[C@@]1(CCC(=CC2)C)C2=C(C(=O)O1)C"),  # Original pattern
        Chem.MolFromSmarts("[C@@]1(CCC(=CC2)C)C(=C)(C(=O)O1)C2"),  # Alternative pattern
        Chem.MolFromSmarts("[C@@]1(CCC(=CC2)C)C2=C(C1=O)C")  # Alternative pattern without alcohol
    ]
    has_germacrane_skeleton = any(mol.HasSubstructMatch(pattern) for pattern in germacrane_patterns)
    if not has_germacrane_skeleton:
        return False, "Does not contain germacrane skeleton"
    
    # Look for lactone ring (-O-C(=O)-)
    lactone_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    has_lactone_ring = mol.HasSubstructMatch(lactone_pattern)
    if not has_lactone_ring:
        return False, "No lactone ring found"
    
    # Check sesquiterpene nature (15 carbon atoms)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 15:
        return False, "Not a sesquiterpene (does not have 15 carbon atoms)"
    
    # Check for common functional groups in germacranolides
    has_alcohol = mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2H]"))
    has_ester = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=[OX1])[OX2]"))
    has_exocyclic_double_bond = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3]=[CX3]"))
    has_ketone = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3]=O"))
    
    if not (has_alcohol or has_ester or has_exocyclic_double_bond or has_ketone):
        return False, "Does not contain common germacranolide functional groups"
    
    return True, "Contains germacrane skeleton with lactone ring and sesquiterpene nature"