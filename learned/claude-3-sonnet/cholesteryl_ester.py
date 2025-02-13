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
        
    # Look for basic steroid core (four fused rings, more permissive pattern)
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for characteristic cholesterol features including:
    # - Double bond at C5-C6 position
    # - Branched aliphatic tail
    # - Connection point for ester at C3
    cholesterol_features = Chem.MolFromSmarts("[#6]-[#6]-[#6]-[#6]-1-[#6]-[#6]-[#6]~2-[#6]~[#6]-[#6]=[#6]-[#6]-[#6]-2")
    if not mol.HasSubstructMatch(cholesterol_features):
        return False, "Missing characteristic cholesterol features"

    # Check for iso-octyl tail characteristic of cholesterol
    tail_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][CH](C)[CH2][CH](C)C")
    if not mol.HasSubstructMatch(tail_pattern):
        return False, "Missing characteristic cholesterol side chain"

    # Verify oxygen count (2 for ester group)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:
        return False, f"Found {oxygen_count} oxygens, should be exactly 2 for cholesteryl ester"

    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450:  # Lowered threshold slightly
        return False, "Molecular weight too low for cholesteryl ester"

    # Count carbons - cholesteryl esters typically have >27 carbons (cholesterol) + fatty acid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 29:
        return False, "Too few carbons for cholesteryl ester"

    # Check for fatty acid chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Acyl chain too short for typical cholesteryl ester"

    return True, "Contains cholesterol core with ester-linked fatty acid"