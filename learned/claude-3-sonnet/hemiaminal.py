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

    # Basic hemiaminal pattern: Carbon with both OH and N attached
    # [CX4] - sp3 carbon
    # -[OX2H1] - hydroxyl group
    # -[#7] - any nitrogen
    basic_pattern = Chem.MolFromSmarts('[CX4](-[OX2H1])(-[#7])')
    
    # Alternative pattern for cyclic structures
    cyclic_pattern = Chem.MolFromSmarts('[CX4;R](-[OX2H1])(-[#7;R])')
    
    # Pattern for bridged structures often found in natural products
    bridged_pattern = Chem.MolFromSmarts('[CX4](-[OX2H1])(-[#7X3])')

    matches = []
    
    # Check all patterns
    for pattern in [basic_pattern, cyclic_pattern, bridged_pattern]:
        if pattern is not None:
            matches.extend(mol.GetSubstructMatches(pattern))

    if not matches:
        return False, "No hemiaminal group found"

    # Exclude false positives
    exclude_patterns = [
        '[NX3](=[OX1])',  # Nitro group
        '[CX3](=O)[OX2H1]',  # Carboxylic acid
        '[NX3]-[OX2]'  # N-O bond (hydroxylamine)
    ]
    
    # Convert matches to set of unique carbons (the first atom in each match)
    carbon_matches = set(match[0] for match in matches)
    
    # Remove carbons that are part of excluded patterns
    for pattern in exclude_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            exclude_matches = mol.GetSubstructMatches(pat)
            for match in exclude_matches:
                carbon_matches.discard(match[0])

    if not carbon_matches:
        return False, "Found potential matches but they are excluded patterns"

    # Additional validation
    valid_carbons = set()
    for carbon_idx in carbon_matches:
        carbon = mol.GetAtomWithIdx(carbon_idx)
        # Check if carbon has both O and N neighbors
        neighbors = carbon.GetNeighbors()
        has_oh = False
        has_n = False
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 8:
                # Check if it's a hydroxyl
                if neighbor.GetTotalNumHs() > 0:
                    has_oh = True
            elif neighbor.GetAtomicNum() == 7:
                has_n = True
        if has_oh and has_n:
            valid_carbons.add(carbon_idx)

    if not valid_carbons:
        return False, "No valid hemiaminal groups found after structural validation"

    reason = f"Found {len(valid_carbons)} valid hemiaminal group(s)"
    return True, reason