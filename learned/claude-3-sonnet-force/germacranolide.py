"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: CHEBI:59245 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import NumAromaticRings, MolWt
from rdkit.Chem.Lipinski import NumRotatableBonds

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
    
    # Check for sesquiterpene (15 carbon atoms)
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if carbon_count != 15:
        return False, "Not a sesquiterpene (does not have 15 carbon atoms)"
    
    # Check for lactone ring
    lactone_pattern = Chem.MolFromSmarts("[C!R]1C(=O)OC1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"
    
    # Check for germacrane skeleton
    germacrane_patterns = [
        Chem.MolFromSmarts("[C]1(C)[C@@](CCC[C@H]2[C@@H]1CC[C@@H]2C)([H])C"),  # Pattern 1
        Chem.MolFromSmarts("[C]1(C)[C@](CCC[C@H]2[C@@H]1CC[C@@H]2C)([H])C"),  # Pattern 2 (different stereochemistry)
        Chem.MolFromSmarts("[C]1(C)[C@@H](CCC[C@H]2[C@@H]1CC[C@]2(C)C)C")  # Pattern 3 (different substituent)
    ]
    
    for pattern in germacrane_patterns:
        if mol.HasSubstructMatch(pattern):
            break
    else:
        return False, "Germacrane skeleton not found"
    
    # Check for common functional groups
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    exocyclic_double_bond_pattern = Chem.MolFromSmarts("[CX3]=C")
    ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX3]")
    
    has_alcohol = mol.HasSubstructMatch(alcohol_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_exocyclic_double_bond = mol.HasSubstructMatch(exocyclic_double_bond_pattern)
    has_ketone = mol.HasSubstructMatch(ketone_pattern)
    
    # Check for molecular weight and rotatable bonds
    mol_wt = MolWt(mol)
    num_rotatable_bonds = NumRotatableBonds(mol)
    
    if mol_wt < 200 or mol_wt > 450:
        return False, "Molecular weight outside the typical range for germacranolides"
    
    if num_rotatable_bonds < 3 or num_rotatable_bonds > 10:
        return False, "Number of rotatable bonds outside the typical range for germacranolides"
    
    # If all conditions are met, classify as germacranolide
    return True, "Molecule contains the germacrane skeleton, a lactone ring, and common functional groups of germacranolides"