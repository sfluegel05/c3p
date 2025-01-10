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

    # Core hemiaminal pattern: carbon with both OH and N attached
    # Excludes cases where N is part of amide/peptide or OH is part of carboxylic acid
    # [CX4] ensures sp3 carbon
    # !$([#6](=O)[OH1]) excludes carboxylic acids
    # !$([#7]C=O) excludes amides
    hemiaminal_pattern = Chem.MolFromSmarts('''
        [CX4;!$(C(=O));!$(C([#7]C=O))]-[OH1]
        AND
        [CX4;!$(C(=O));!$(C([#7]C=O))]-[#7;!$(NC=O)]
    ''')

    # Additional pattern for cyclic hemiaminals
    cyclic_hemiaminal_pattern = Chem.MolFromSmarts('''
        [CX4;!$(C(=O));!$(C([#7]C=O))](-[OH1])-[#7;R;!$(NC=O)]
    ''')

    # Exclusion patterns
    exclude_patterns = [
        Chem.MolFromSmarts('[NX3]-[OX2]'),  # N-hydroxy
        Chem.MolFromSmarts('[CX3](=O)[OX2]'),  # carboxyl/ester
        Chem.MolFromSmarts('[CX3](=O)[NX3]'),  # amide
        Chem.MolFromSmarts('[C]-[N+](-[O-])')  # N-oxide
    ]

    # Get matches
    matches = []
    if hemiaminal_pattern:
        matches.extend(mol.GetSubstructMatches(hemiaminal_pattern))
    if cyclic_hemiaminal_pattern:
        matches.extend(mol.GetSubstructMatches(cyclic_hemiaminal_pattern))

    if not matches:
        return False, "No hemiaminal group found"

    # Filter out matches that overlap with exclusion patterns
    exclude_atoms = set()
    for pattern in exclude_patterns:
        if pattern:
            for match in mol.GetSubstructMatches(pattern):
                exclude_atoms.update(match)

    # Keep only valid matches
    valid_matches = []
    for match in matches:
        if not any(atom in exclude_atoms for atom in match):
            valid_matches.append(match)

    if not valid_matches:
        return False, "Found potential matches but they are excluded patterns"

    # Additional validation: verify OH and N are on same carbon
    confirmed_matches = []
    for match in valid_matches:
        carbon = mol.GetAtomWithIdx(match[0])
        if carbon.GetHybridization() == Chem.HybridizationType.SP3:
            # Check neighbors for both OH and N
            has_oh = False
            has_n = False
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:
                    has_oh = True
                elif neighbor.GetAtomicNum() == 7 and not any(n.GetAtomicNum() == 8 and n.GetBonds()[0].GetBondType() == Chem.BondType.DOUBLE 
                                                            for n in neighbor.GetNeighbors()):
                    has_n = True
            if has_oh and has_n:
                confirmed_matches.append(match)

    if not confirmed_matches:
        return False, "No valid hemiaminal groups found after structural validation"

    reason = f"Found {len(confirmed_matches)} valid hemiaminal group(s)"
    return True, reason