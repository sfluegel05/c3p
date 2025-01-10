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

    # Define patterns for amino groups and groups to exclude
    patterns = {
        # Primary amines (aliphatic and aromatic)
        'primary_amine': '[NX3H2][CX4,cX3]',
        
        # Secondary amines (not in rings)
        'secondary_amine': '[NX3H1;!R]([CX4,cX3])[CX4,cX3]',
        
        # Secondary amines (in rings, but not aromatic)
        'cyclic_amine': '[NX3H1;R;!a]([CX4,cX3])[CX4,cX3]',
        
        # Substituted amines (including N-hydroxy)
        'substituted_amine': '[NX3;!$(NC=O);!$(N=*);!$(N#*)]([CX4,cX3])',
        
        # Protonated amines
        'protonated_amine': '[NX4H3+,NX4H2+,NX4H+]',
        
        # Amino groups in triazines (if not part of amide)
        'triazine_amine': '[NX3;!$(NC=O);!$(N=*)]([c])',
        
        # Patterns to exclude
        'exclude_amide': '[NX3][CX3](=[OX1])',
        'exclude_nitro': '[$([NX3](=O)=O),$([NX3+](=O)[O-])]',
        'exclude_azo': '[NX2]=[NX2]',
        'exclude_nitrile': '[NX1]#[CX2]',
        'exclude_aromatic': '[nX2,nX3;a]', # aromatic N in 5 or 6 membered rings
        'exclude_imine': '[NX2]=[CX3]',
        'exclude_hydrazine': '[NX3H2]-[NX3H2]'
    }

    # Convert patterns to RDKit molecules
    patterns = {name: Chem.MolFromSmarts(pattern) for name, pattern in patterns.items()}

    # Get all matches for amino groups
    amino_matches = set()
    
    for pattern_name in ['primary_amine', 'secondary_amine', 'cyclic_amine', 
                        'substituted_amine', 'protonated_amine', 'triazine_amine']:
        if patterns[pattern_name] is not None:
            matches = mol.GetSubstructMatches(patterns[pattern_name])
            for match in matches:
                amino_matches.add(match[0])  # Add nitrogen atom index

    # Remove matches that are part of excluded groups
    for exclude_pattern in ['exclude_amide', 'exclude_nitro', 'exclude_azo', 
                          'exclude_nitrile', 'exclude_aromatic', 'exclude_imine',
                          'exclude_hydrazine']:
        if patterns[exclude_pattern] is not None:
            exclude_matches = mol.GetSubstructMatches(patterns[exclude_pattern])
            for match in exclude_matches:
                if match[0] in amino_matches:
                    amino_matches.remove(match[0])

    total_amino_groups = len(amino_matches)

    if total_amino_groups >= 2:
        return True, f"Contains {total_amino_groups} amino groups"
    
    return False, f"Contains only {total_amino_groups} amino groups, minimum 2 required"