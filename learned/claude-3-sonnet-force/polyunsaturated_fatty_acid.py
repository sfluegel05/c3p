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
    
    # Check for long aliphatic chain with multiple double bonds
    pufa_pattern = Chem.MolFromSmarts("C(=O)OCCCC(/C=C/C)C(/C=C/C)C")
    if not mol.HasSubstructMatch(pufa_pattern):
        return False, "Does not match the pattern of a polyunsaturated fatty acid"
    
    # Count number of double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    num_double_bonds = len(double_bond_matches)
    if num_double_bonds < 2:
        return False, "Less than two double bonds found"
    
    # Check for cis/trans configuration of double bonds
    cis_trans_pattern = Chem.MolFromSmarts("/C=C/C")
    cis_trans_matches = mol.GetSubstructMatches(cis_trans_pattern)
    if len(cis_trans_matches) == 0:
        return False, "No cis or trans double bonds found"
    
    # Additional checks or constraints can be added here, if needed
    
    return True, f"Contains {num_double_bonds} double bonds in a long aliphatic chain with a carboxylic acid group"