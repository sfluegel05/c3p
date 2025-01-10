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
    
    # Multiple SMARTS patterns for different types of indole cores
    indole_patterns = [
        'c12ccccc1[nH]cc2',           # Basic indole
        'c12ccccc1[nX3]cc2',          # N-substituted indole
        'c12ccccc1nc(C)c2',           # Alternative indole pattern
        'c12ccccc1[nX3]c(C)c2',       # N-substituted with C-substitution
        'c12ccccc1[nX3]C=C2',         # Modified indole double bond
        'c12ccccc1[nX3]CC2',          # Reduced indole
        'c12ccccc1[nX3][C,c]2',       # Bridged indole
        'c12ccccc1n([*])cc2',         # Any N-substituted indole
        'c12ccccc1n([*])c([*])c2',    # Heavily substituted indole
        'C1=CC2=C(C=C1)N=C([*])C2',   # Modified indole with imine
        'C1=CC2=C(C=C1)NC=C2',        # Basic indole alternative
        'C1=CC2=C(C=C1)N([*])C=C2'    # N-substituted alternative
    ]
    
    has_indole = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) 
                    for pattern in indole_patterns)
    
    if not has_indole:
        return False, "No indole core structure found"

    # Check for nitrogen content (alkaloid requirement)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count == 0:
        return False, "No nitrogen atoms found - not an alkaloid"
    
    # Structural checks typical for indole alkaloids
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, "Too few rings for an indole alkaloid"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for an indole alkaloid"
    
    # Carbon count check
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, "Too few carbons for an indole alkaloid"
    
    # Check nitrogen connectivity
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