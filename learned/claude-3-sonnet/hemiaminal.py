"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:60324 hemiaminal
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule contains a hemiaminal group based on its SMILES string.
    A hemiaminal has both an amino group and a hydroxy group attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a hemiaminal group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # Collection of SMARTS patterns for different hemiaminal types
    patterns = [
        # Basic hemiaminal pattern - more relaxed version
        '[CX4](-[OX2])(-[NX3])',
        
        # Cyclic hemiaminal patterns
        '[CX4;R](-[OX2])(-[NX3;R])',
        '[CX4;R](-[OX2])(-[NX3])',
        
        # Pattern for charged species
        '[CX4](-[OX2])(-[NX4+])',
        
        # Bridged bicyclic patterns common in natural products
        '[CX4;R2](-[OX2])(-[NX3;R])',
        
        # Pattern for fused ring systems
        '[CX4;R2](-[OX2;R0])(-[NX3;R])',
        
        # Pattern catching more complex cases
        '[CX4](-[OX2])(-[NX3,NX4+])',
    ]

    matches = set()
    for pattern in patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None:
            these_matches = mol.GetSubstructMatches(pat)
            matches.update(match[0] for match in these_matches)  # Add carbon indices

    if not matches:
        return False, "No hemiaminal group found"

    # Patterns to exclude false positives
    exclude_patterns = [
        '[NX3](=[OX1])',  # Nitro group
        '[CX3](=O)[OX2H1]',  # Carboxylic acid
        '[NX3]-[OX2]',  # N-O bond
        '[CX3](=O)-[OX2]',  # Ester
        '[CX3]=N',  # Imine
        '[CX3](=[OX1])[NX3]'  # Amide
    ]
    
    # Remove false positives
    for pattern in exclude_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            exclude_matches = mol.GetSubstructMatches(pat)
            for match in exclude_matches:
                matches.discard(match[0])

    if not matches:
        return False, "Found potential matches but they are excluded patterns"

    # Additional validation
    valid_carbons = set()
    for carbon_idx in matches:
        carbon = mol.GetAtomWithIdx(carbon_idx)
        neighbors = carbon.GetNeighbors()
        
        # Check neighbor atoms
        n_count = 0
        o_count = 0
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 7:  # N
                if neighbor.GetFormalCharge() >= 0:  # Allow neutral or positive N
                    n_count += 1
            elif neighbor.GetAtomicNum() == 8:  # O
                if neighbor.GetFormalCharge() == 0:  # Neutral O only
                    o_count += 1
        
        # Must have exactly one N and one O
        if n_count == 1 and o_count == 1:
            # Check hybridization of carbon
            if carbon.GetHybridization() == Chem.HybridizationType.SP3:
                valid_carbons.add(carbon_idx)

    if not valid_carbons:
        return False, "No valid hemiaminal groups found after structural validation"

    reason = f"Found {len(valid_carbons)} valid hemiaminal group(s)"
    return True, reason