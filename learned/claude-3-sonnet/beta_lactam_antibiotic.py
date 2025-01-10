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

    # More general beta-lactam ring pattern
    # Allows for different ring fusions and substitution patterns
    beta_lactam_pattern = Chem.MolFromSmarts("[N]1[C]([C][C]1)=O")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Common structural features of beta-lactam antibiotics
    patterns = {
        # Penicillin-like structure (thiazolidine ring fused to beta-lactam)
        'penam': Chem.MolFromSmarts("[N]1[C]2[C][S][C][C]2[C]1=O"),
        
        # Cephalosporin-like structure
        'cephem': Chem.MolFromSmarts("[N]1[C]2[C][S][C][C]2[C]1=O"),
        
        # Carbapenem-like structure
        'carbapenem': Chem.MolFromSmarts("[N]1[C]2[C][C][C]2[C]1=O"),
        
        # Monobactam-like structure with sulfonic acid
        'monobactam': Chem.MolFromSmarts("[N]1[C]([C][C]1)=O.[S](=O)(=O)[O]"),
        
        # Oxacephem-like structure
        'oxacephem': Chem.MolFromSmarts("[N]1[C]2[C][O][C][C]2[C]1=O"),
        
        # Carboxylic acid group (can be deprotonated)
        'carboxyl': Chem.MolFromSmarts("C(=O)[OH,O-]"),
        
        # Common amide side chain
        'amide': Chem.MolFromSmarts("NC(=O)")
    }

    # Check for required structural features
    matches = {name: mol.HasSubstructMatch(pattern) for name, pattern in patterns.items()}
    
    # Must have carboxylic acid group
    if not matches['carboxyl']:
        return False, "No carboxylic acid group found"

    # Must have at least one of the common beta-lactam ring systems
    ring_types = ['penam', 'cephem', 'carbapenem', 'monobactam', 'oxacephem']
    if not any(matches[ring_type] for ring_type in ring_types):
        if not matches['amide']:  # If no typical ring system, at least require amide
            return False, "Missing characteristic beta-lactam antibiotic structural features"

    # Basic molecular properties checks
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 10:
        return False, "Molecule too small to be a beta-lactam antibiotic"

    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for beta-lactam antibiotic"

    # Count nitrogen atoms (beta-lactams typically have multiple nitrogens)
    n_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7])
    if n_count < 1:
        return False, "Insufficient nitrogen atoms for beta-lactam antibiotic"

    # Look for characteristic ring fusion patterns
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "No rings found"

    return True, "Contains beta-lactam ring and characteristic features of beta-lactam antibiotics"