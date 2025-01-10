"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:26144 polyamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is an organic compound containing two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is organic (contains carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not an organic compound"

    # Define patterns for different types of amino groups
    patterns = {
        # Primary amines (aliphatic and aromatic)
        'primary_amine': '[NX3H2]',
        
        # Secondary amines (aliphatic and aromatic)
        'secondary_amine': '[NX3H1]([#6,#1])[#6]',
        
        # Tertiary amines (aliphatic and aromatic)
        'tertiary_amine': '[NX3H0]([#6])[#6]',
        
        # Protonated amines
        'protonated_amine': '[NX4H3+,NX4H2+,NX4H+]',
        
        # Triazine amino groups
        'triazine_amine': '[NX3H2,NX3H1,NX3H0]([#6])[#6]',
        
        # Patterns to exclude
        'amide': '[NX3][CX3](=[OX1])',
        'nitro': '[$([NX3](=O)=O),$([NX3+](=O)[O-])]',
        'azo': '[NX2]=[NX2]',
        'nitrile': '[NX1]#[CX2]',
    }

    # Convert patterns to RDKit molecules
    patterns = {name: Chem.MolFromSmarts(pattern) for name, pattern in patterns.items()}

    # Count different types of amino groups
    counts = {}
    for name, pattern in patterns.items():
        matches = mol.GetSubstructMatches(pattern)
        counts[name] = len(matches)

    # Get all nitrogen atoms
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    
    # Count amino groups (excluding those in excluded patterns)
    amino_nitrogens = set()
    
    for n in nitrogens:
        # Skip if nitrogen is part of excluded groups
        is_excluded = False
        for exclude_pattern in ['amide', 'nitro', 'azo', 'nitrile']:
            matches = mol.GetSubstructMatches(patterns[exclude_pattern])
            for match in matches:
                if n.GetIdx() in match:
                    is_excluded = True
                    break
        if is_excluded:
            continue
            
        # Check if nitrogen has hydrogens or is part of valid amino group
        if (n.GetTotalNumHs() > 0 or  # Has hydrogens
            n.GetFormalCharge() == 1 or  # Protonated
            any(n.GetIdx() in match for match in mol.GetSubstructMatches(patterns['triazine_amine']))  # Triazine amino
           ):
            amino_nitrogens.add(n.GetIdx())

    total_amino_groups = len(amino_nitrogens)

    if total_amino_groups >= 2:
        return True, f"Contains {total_amino_groups} amino groups"
    
    return False, f"Contains only {total_amino_groups} amino groups, minimum 2 required"