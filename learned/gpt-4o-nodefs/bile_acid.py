"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid is generally characterized by a steroid core,
    with various hydroxyl or keto groups, and a carboxylic acid group in the side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible steroid core pattern using SMARTS for detecting three fused rings with specific bonding
    steroid_core_pattern = Chem.MolFromSmarts("[#6]1([#6])[#6]2[#6]3[#6][#6][#6][#6]3[#6][#6]2[#6][#6]1")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"

    # Check for hydroxyl groups, common in bile acids
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:  # Generally, bile acids may have at least two hydroxyl groups
        return False, "Insufficient hydroxyl groups found"

    # Check for keto groups, less common but present in some bile acids
    keto_pattern = Chem.MolFromSmarts("C(=O)")
    keto_matches = mol.GetSubstructMatches(keto_pattern)

    # Look for a carboxylic acid tail, typical of bile acids
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Contains a steroid core with appropriate substituents including hydroxyl and carboxylic acid groups"