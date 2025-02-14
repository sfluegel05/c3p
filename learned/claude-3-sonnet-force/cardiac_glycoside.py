"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: CHEBI:28807 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    A cardiac glycoside is a steroid lactone containing sugar residues that acts on the contractile force of cardiac muscles.

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
    
    # Look for steroid backbone pattern with lactone ring
    steroid_pattern = Chem.MolFromSmarts("[C@@]12[C@H](CC[C@@]3([C@@H]1CC[C@@]4([C@]3(CC[C@@]4(C(=O)O)C)C)C)C)CC[C@]2(C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone with lactone ring found"
    
    # Look for sugar residues
    sugar_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][OX2][CX4][OX2]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar residues found"
    
    # Count rotatable bonds to verify presence of sugar chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Not enough rotatable bonds for sugar chains"
    
    # Check molecular weight - cardiac glycosides typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for cardiac glycoside"
    
    # Count oxygen atoms - should have several from sugar residues
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 10:
        return False, "Too few oxygen atoms for cardiac glycoside"
    
    return True, "Contains steroid backbone with lactone ring and sugar residues"