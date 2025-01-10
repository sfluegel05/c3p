"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: indole alkaloid
Definition: An alkaloid containing an indole skeleton.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for indole core structure
    indole_pattern = Chem.MolFromSmarts('c12ccccc1[nH]cc2')  # Basic indole skeleton
    indole_pattern_2 = Chem.MolFromSmarts('c12ccccc1nc(C)c2')  # Alternative indole pattern
    indole_pattern_3 = Chem.MolFromSmarts('c12ccccc1[nX3]cc2')  # N-substituted indole
    
    has_indole = mol.HasSubstructMatch(indole_pattern) or \
                 mol.HasSubstructMatch(indole_pattern_2) or \
                 mol.HasSubstructMatch(indole_pattern_3)
    
    if not has_indole:
        return False, "No indole core structure found"

    # Count nitrogens to confirm it's an alkaloid
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count == 0:
        return False, "No nitrogen atoms found - not an alkaloid"
        
    # Check for typical characteristics of alkaloids
    # Most indole alkaloids have multiple rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, "Too few rings for an indole alkaloid"
        
    # Most natural indole alkaloids have significant molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for an indole alkaloid"
        
    # Check for presence of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, "Too few carbons for an indole alkaloid"
        
    # Additional check for connectivity
    # Most indole alkaloids have the nitrogen as part of a larger system
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    for n_atom in n_atoms:
        if len(n_atom.GetNeighbors()) >= 2:
            return True, "Contains indole core and has characteristics of an alkaloid"
            
    return False, "Nitrogen not properly incorporated into ring system"

__metadata__ = {
    'chemical_class': {
        'name': 'indole alkaloid',
        'definition': 'An alkaloid containing an indole skeleton.',
        'parents': ['alkaloid', 'indole-containing compound']
    }
}