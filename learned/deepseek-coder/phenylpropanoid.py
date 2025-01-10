"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: CHEBI:26195 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is an organic aromatic compound with a structure based on a phenylpropane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one benzene ring
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"

    # More flexible phenylpropane pattern matching
    # Allows for variations in the propyl chain and different substitutions
    phenylpropane_patterns = [
        # Standard patterns
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4]"),  # Standard phenylpropane
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4]=[CX4]"),      # Double bond in chain
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4]"),       # Shorter chain
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][CX4]"),  # Longer chain
        
        # More flexible patterns with substitutions
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][*]"),  # Substituted chain
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4]=[CX4][*]"),      # Substituted double bond
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][*]"),       # Substituted shorter chain
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][CX4][*]"),  # Substituted longer chain
        
        # Patterns with different connections
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][c]"),  # Connected to another aromatic
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][O]"),  # Connected to oxygen
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][N]"),  # Connected to nitrogen
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][S]"),  # Connected to sulfur
    ]

    has_phenylpropane = any(mol.HasSubstructMatch(pattern) for pattern in phenylpropane_patterns)
    if not has_phenylpropane:
        return False, "No phenylpropane-like skeleton found"

    # Check molecular weight - phenylpropanoids typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for phenylpropanoid"

    # Additional check for common phenylpropanoid features
    # Look for common functional groups in phenylpropanoids
    common_groups = [
        Chem.MolFromSmarts("[OX2H]"),  # Hydroxyl
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),  # Carboxylic acid
        Chem.MolFromSmarts("[CX3](=O)[OX2]"),  # Ester
        Chem.MolFromSmarts("[CX3](=O)[NX3]"),  # Amide
    ]
    
    has_common_groups = any(mol.HasSubstructMatch(group) for group in common_groups)
    if not has_common_groups:
        return False, "No common phenylpropanoid functional groups found"

    # If all checks pass, classify as phenylpropanoid
    return True, "Contains phenylpropane-like skeleton, common functional groups, and meets molecular weight requirements"