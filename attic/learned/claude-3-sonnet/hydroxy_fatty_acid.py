"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: hydroxy fatty acid
A fatty acid carrying one or more hydroxy substituents
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for hydroxyl groups (excluding the one in COOH)
    # First find all OH groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Then find COOH groups to subtract
    cooh_matches = mol.GetSubstructMatches(carboxylic_pattern)
    
    # Count hydroxyls that are not part of COOH
    non_acid_hydroxyls = len(hydroxyl_matches) - len(cooh_matches)
    
    if non_acid_hydroxyls < 1:
        return False, "No hydroxyl groups found (excluding carboxylic acid)"
    
    # Check carbon chain length (should be at least 3 carbons for fatty acid)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, "Carbon chain too short to be a fatty acid"
    
    # Check for continuous carbon chain
    alkyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "No continuous carbon chain found"
    
    # Additional check for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 60:  # Minimum weight for a hydroxy fatty acid
        return False, "Molecular weight too low for hydroxy fatty acid"
        
    # Count rotatable bonds to verify chain nature
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 1:
        return False, "Too rigid to be a fatty acid"
    
    return True, f"Contains carboxylic acid group and {non_acid_hydroxyls} additional hydroxyl group(s) on carbon chain"