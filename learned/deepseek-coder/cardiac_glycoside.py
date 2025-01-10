"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: CHEBI:36218 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    A cardiac glycoside is a steroid lactone containing sugar residues that act on the contractile force of the cardiac muscles.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid nucleus (cyclopentanoperhydrophenanthrene ring system)
    # More general pattern to match different configurations
    steroid_pattern = Chem.MolFromSmarts("[C]12[C][C][C]3[C]([C]1[C][C][C]2[C])[C]4[C][C][C]([C]4(C)[C]3)O")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"

    # Check for lactone ring (butenolide or pentenolide)
    # More specific pattern for lactone rings
    lactone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4][CX4]=[CX3]")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Check for sugar residues (glycosidic bonds)
    sugar_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 1:
        return False, "No sugar residues found"

    # Check molecular weight - cardiac glycosides typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for cardiac glycoside"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for cardiac glycoside"
    if o_count < 6:
        return False, "Too few oxygens for cardiac glycoside"

    return True, "Contains steroid nucleus, lactone ring, and sugar residues"