"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: cholesteryl ester
A sterol ester obtained by formal condensation of the carboxy group of any 
carboxylic acid with the 3-hydroxy group of cholesterol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
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

    # Look for steroid core (four fused rings) - simplified pattern
    steroid_core = Chem.MolFromSmarts("C1CC2CCC3C(C2)CCC4CCCC34C1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for double bond characteristic of cholesterol (typically at C5-C6)
    double_bond = Chem.MolFromSmarts("C=CC1CCC2C(C1)CCC")
    if not mol.HasSubstructMatch(double_bond):
        return False, "Missing characteristic cholesterol double bond"

    # Check for branched side chain characteristic of cholesterol
    # Pattern matches iso-octyl tail: -CH2-CH2-CH2-CH(CH3)-CH2-CH(CH3)-CH3
    side_chain = Chem.MolFromSmarts("CCC(C)CCC(C)C")
    if not mol.HasSubstructMatch(side_chain):
        return False, "Missing characteristic cholesterol side chain"

    # Verify oxygen count (exactly 2 for ester group)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:
        return False, f"Found {oxygen_count} oxygens, should be exactly 2 for cholesteryl ester"

    # Count carbons - cholesteryl esters typically have >27 carbons 
    # (27 from cholesterol + at least 2 from acyl group)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 29:
        return False, "Too few carbons for cholesteryl ester"

    # Check for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:  # Lowered threshold to catch smaller esters
        return False, "Molecular weight too low for cholesteryl ester"

    # Check for fatty acid chain length via rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Acyl chain too short for typical cholesteryl ester"

    return True, "Contains cholesterol core with ester-linked fatty acid chain"