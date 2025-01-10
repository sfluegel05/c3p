"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: CHEBI:35757 tricarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is an oxoacid containing three carboxy groups (-C(=O)OH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxylic acid group pattern
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    
    # Find all matches of the carboxylic acid group
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    
    # Check if there are exactly 3 carboxylic acid groups
    if len(carboxy_matches) != 3:
        return False, f"Found {len(carboxy_matches)} carboxylic acid groups, need exactly 3"

    # Ensure the molecule is not a peptide or complex structure
    # by checking the number of rotatable bonds and molecular weight
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Tricarboxylic acids typically have a moderate number of rotatable bonds
    # and a molecular weight within a certain range
    if n_rotatable > 20 or mol_wt > 1000:
        return False, "Molecule is too complex or too large to be a simple tricarboxylic acid"

    # Check that the carboxylic acid groups are not part of a peptide
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    if len(peptide_matches) > 0:
        return False, "Molecule contains peptide bonds, not a simple tricarboxylic acid"

    return True, "Contains exactly 3 carboxylic acid groups and is not a complex structure"