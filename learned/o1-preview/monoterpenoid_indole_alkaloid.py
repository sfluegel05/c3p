"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    A monoterpenoid indole alkaloid contains an indole moiety derived from L-tryptophan and
    a monoterpene-derived moiety (C10 unit), usually linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for indole moiety (allowing substitutions at nitrogen)
    indole_pattern = Chem.MolFromSmarts('c1ccc2c(c1)[nH]c[cH]2')  # Indole core
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole moiety found"
        
    # Check for monoterpene-derived moiety (approximate)
    # Monoterpenes are C10 units built from isoprene units
    # We can attempt to find fragments with 10 carbons (C10)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 20:
        return False, f"Too few carbons ({num_carbons}) to be a monoterpenoid indole alkaloid"
    
    # Check for isoprene units (C5H8) repeated twice
    isoprene_pattern = Chem.MolFromSmarts('C=C(C)C')
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, "Less than two isoprene units found"
        
    # Check for connectivity between indole and monoterpene units
    # This is a complex task and may require detailed analysis
    # For simplicity, we'll assume presence of both moieties is sufficient
    
    return True, "Contains indole moiety and monoterpenoid unit"

__metadata__ = {
    'chemical_class': {
        'name': 'monoterpenoid indole alkaloid',
        'definition': 'A terpenoid indole alkaloid which is biosynthesised from L-tryptophan and diisoprenoid (usually secologanin) building blocks.'
    },
    'config': {
        # Configuration parameters if any
    }
}