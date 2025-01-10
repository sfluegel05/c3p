"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: beta-lactam antibiotics
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    Beta-lactam antibiotics contain a beta-lactam ring and are used as antibiotics.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove any salt counterions to analyze just the main molecule
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(fragments) > 1:
        # Take the largest fragment as the main molecule
        mol = max(fragments, key=lambda x: x.GetNumAtoms())

    # Essential beta-lactam ring pattern (4-membered ring with N and C=O)
    beta_lactam_pattern = Chem.MolFromSmarts("[NR1]1[CR2][CR2][CR2]1(=O)")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Common structural patterns in beta-lactam antibiotics
    patterns = {
        # Penicillin-like structure (thiazolidine ring fused to beta-lactam)
        'penam': Chem.MolFromSmarts("[NR1]1[C@@H]2[C@H]([C]1=O)[SC](C)(C)[CH]2"),
        
        # Cephalosporin-like structure
        'cephem': Chem.MolFromSmarts("[NR1]1[C]2[CH][SC]=C2[C]1=O"),
        
        # Carbapenem-like structure
        'carbapenem': Chem.MolFromSmarts("[NR1]1[C]2[CH][CH][C@H]2[C]1=O"),
        
        # Monobactam-like structure
        'monobactam': Chem.MolFromSmarts("[NR1]1[CH][CH][C]1=O"),
        
        # Oxacephem-like structure
        'oxacephem': Chem.MolFromSmarts("[NR1]1[C]2[CH][OC][C]2[C]1=O"),
        
        # Carboxylic acid or derivative group
        'carboxyl': Chem.MolFromSmarts("[$([C](=O)[OH]),$([C](=O)[O-]),$([C](=O)O[CH2,CH3])]"),
        
        # Amide side chain
        'amide': Chem.MolFromSmarts("NC(=O)")
    }

    # Check for structural features
    matches = {name: mol.HasSubstructMatch(pattern) for name, pattern in patterns.items()}
    
    # Must have carboxylic acid or derivative group
    if not matches['carboxyl']:
        return False, "No carboxylic acid or derivative group found"

    # Must have at least one of the common beta-lactam ring systems
    ring_types = ['penam', 'cephem', 'carbapenem', 'monobactam', 'oxacephem']
    has_characteristic_ring = any(matches[ring_type] for ring_type in ring_types)
    
    # Count rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    # Count atoms
    num_atoms = mol.GetNumAtoms()
    
    # Count nitrogens (excluding charges)
    n_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7])
    
    # Special cases for known beta-lactam cores
    if num_atoms < 7:
        return False, "Molecule too small to be a beta-lactam antibiotic"
        
    if has_characteristic_ring:
        if num_rings >= 2 and n_count >= 1:
            return True, "Contains characteristic beta-lactam antibiotic ring system"
    else:
        # For non-standard structures, require more evidence
        if num_rings >= 2 and n_count >= 2 and matches['amide']:
            return True, "Contains beta-lactam ring with antibiotic-like features"
            
    # For simple beta-lactam cores
    if matches['monobactam'] and matches['carboxyl']:
        return True, "Contains basic beta-lactam antibiotic structure"

    return False, "Missing required structural features for beta-lactam antibiotic classification"