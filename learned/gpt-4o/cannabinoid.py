"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a phenolic group - common in cannabinoids
    phenol_pattern = Chem.MolFromSmarts('c1cc(O)ccc1')
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenolic group detected - common in cannabinoids"

    # Check for long carbon chains, typically present in cannabinoids
    chain_pattern = Chem.MolFromSmarts('CCCCCCCCCCCC')  # Long C-chain
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long carbon chain detected, typical of cannabinoids"

    # Check for oxygen in heterocyclic ring or as functional group
    oxygen_pattern = Chem.MolFromSmarts('[O]')
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "No heterocyclic oxygen or oxygen-containing group found"

    # These checks confirm essential cannabinoid characteristics
    return True, "Phenolic group, long carbon chain, and oxygen-containing functional groups detected"