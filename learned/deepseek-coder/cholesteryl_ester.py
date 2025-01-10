"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: CHEBI:36235 cholesteryl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is a sterol ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of cholesterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible cholesterol backbone pattern
    # This pattern matches the core structure of cholesterol, allowing for different stereochemistry and hydrogen atoms
    cholesterol_pattern = Chem.MolFromSmarts("[C@@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@H]3CC[C@@H]4[C@@H](C)CCCC(C)C)")
    if not mol.HasSubstructMatch(cholesterol_pattern):
        return False, "No cholesterol backbone found"

    # Look for ester group (-O-C(=O)-) attached to the 3-hydroxy group of cholesterol
    ester_pattern = Chem.MolFromSmarts("[CX4][OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester group found attached to the cholesterol backbone"

    # Check for a fatty acid chain (long carbon chain attached to the ester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "No fatty acid chain found attached to the ester"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - cholesteryl esters typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for cholesteryl ester"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for cholesteryl ester"
    if o_count < 2:
        return False, "Must have at least 2 oxygens (ester group)"

    return True, "Contains cholesterol backbone with an ester-linked fatty acid chain"