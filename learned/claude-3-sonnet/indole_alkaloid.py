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
    
    # Core indole patterns and their variations commonly found in alkaloids
    indole_patterns = [
        # Basic indole cores
        'c12ccccc1[nH]cc2',           # Basic indole
        'c12ccccc1[nX3]cc2',          # N-substituted indole
        'c12ccccc1nc([*])c2',         # Alternative indole
        
        # Complex fused systems
        'c12ccccc1n([*])c([*])([*])c2',  # Highly substituted
        'C1=CC2=C(C=C1)N([*])C([*])([*])C2([*])', # Spiro/bridged
        'C1=CC2=C(C=C1)N([*])C3([*])C2([*])C3',   # Bridged systems
        
        # Ergot-type cores
        'C1=CC2=C(C=C1)N([*])C3=C2C([*])=C([*])C3',
        'C1=CC2=C(C=C1)N([*])C3=C2CCN3',
        
        # Iboga-type cores
        'C1=CC2=C(C=C1)N([*])C3C2([*])CCN3',
        'C1=CC2=C(C=C1)N([*])C3=C2C([*])CN3',
        
        # Modified cores common in alkaloids
        'C1=CC=C2C(=C1)N=C([*])C2([*])',
        'C1=CC=C2C(=C1)N([*])=C([*])C2([*])',
        'C1=CC=C2C(=C1)N([*])C(=O)C2([*])',
        
        # Quaternary nitrogen variations
        'c12ccccc1[n+]([*])cc2',
        'C1=CC2=C(C=C1)N+=C([*])C2([*])',
        
        # Reduced forms
        'C1=CC2=C(C=C1)N([*])CC2([*])',
        'C1=CC2=C(C=C1)N([*])C([*])C2([*])',
        
        # Yohimbine-type cores
        'C1=CC2=C(C=C1)N([*])C3([*])C2([*])N([*])C3',
        'C1=CC2=C(C=C1)N([*])C3([*])C2([*])CC3',
        
        # Additional complex systems
        'C1=CC2=C(C=C1)N3C([*])C2([*])C3([*])',
        'C1=CC2=C(C=C1)N([*])C3=C2N([*])C3',
        'C1=CC2=C(C=C1)N([*])C3=C2C(=O)N3',
        'C1=CC2=C(C=C1)N([*])C(=N([*]))C2'
    ]
    
    has_indole = False
    for pattern in indole_patterns:
        try:
            pattern_mol = Chem.MolFromSmarts(pattern)
            if pattern_mol and mol.HasSubstructMatch(pattern_mol):
                has_indole = True
                break
        except:
            continue
    
    if not has_indole:
        return False, "No indole core structure found"

    # Check for nitrogen content
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count == 0:
        return False, "No nitrogen atoms found - not an alkaloid"
    
    # Structural requirements
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    if ring_count < 2:
        return False, "Too few rings for an indole alkaloid"
    
    # Carbon count check
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, "Too few carbons for an indole alkaloid"
    
    # Check for fused ring systems
    ring_systems = []
    bonds = mol.GetBonds()
    for bond in bonds:
        if bond.IsInRing():
            ring_atoms = set()
            for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
                if atom.IsInRing():
                    ring_atoms.add(atom.GetIdx())
            if len(ring_atoms) > 0:
                ring_systems.append(ring_atoms)
    
    # Look for nitrogen atoms in rings and complex ring systems
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    for n_atom in n_atoms:
        if n_atom.IsInRing():
            return True, "Contains indole core and has characteristics of an alkaloid"
    
    return False, "Structure lacks typical indole alkaloid characteristics"

__metadata__ = {
    'chemical_class': {
        'name': 'indole alkaloid',
        'definition': 'An alkaloid containing an indole skeleton.',
        'parents': ['alkaloid', 'indole-containing compound']
    }
}