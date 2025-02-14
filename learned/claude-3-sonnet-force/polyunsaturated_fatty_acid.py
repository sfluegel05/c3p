"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: CHEBI:36344 polyunsaturated fatty acid
Any fatty acid containing more than one double bond. Acids in this group are reported to have cardioprotective effects; and levels are lowered in chronic fatigue syndrome.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for long aliphatic chain
    aliphatic_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")
    if not mol.HasSubstructMatch(aliphatic_chain_pattern):
        return False, "Aliphatic chain too short"
    
    # Check for multiple double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 2:
        return False, "Less than two double bonds found"
    
    # Check position of double bonds relative to carboxylic acid group
    aliphatic_chain_length = rdMolDescriptors.CalcMolWt(mol) // 14  # approximate chain length
    double_bond_positions = [atom_idx for match in double_bond_matches for atom_idx in match]
    distance_from_acid = min(mol.GetAtomWithIdx(idx).GetShortestPathDistanceToSubmits()[0] for idx in double_bond_positions)
    if distance_from_acid < 3 or distance_from_acid > aliphatic_chain_length - 3:
        return False, "Double bonds too close to carboxylic acid group"
    
    # Additional checks or constraints can be added here, if needed
    
    return True, "Contains multiple double bonds in a long aliphatic chain with a carboxylic acid group"