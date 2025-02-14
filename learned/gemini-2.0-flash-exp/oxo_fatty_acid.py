"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is a fatty acid (long chain carboxylic acid) containing at least one
    additional ketone or aldehyde group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
       return False, "No carboxylic acid group found"
    
    # Check for ketone or aldehyde group
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4]")  # Simplified ketone pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CH1X3](=O)") # Correct aldehyde pattern
    
    has_ketone = mol.HasSubstructMatch(ketone_pattern)
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    
    if not (has_ketone or has_aldehyde):
       return False, "No ketone or aldehyde group found"

    # Check for linear carbon chain (at least 4 carbons)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
       return False, "No long carbon chain found"
       
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
         return False, "Too few carbons for a fatty acid"

    # Count rotatable bonds (relaxed constraints)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 1: # Lower threshold
        return False, "Too few rotatable bonds for a fatty acid"

    # Molecular weight check for fatty acid (relaxed constraint)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100: # Lower minimal MW
      return False, "Molecular weight too low for fatty acid"

    return True, "Contains a carboxylic acid, a ketone or aldehyde, and a long carbon chain"