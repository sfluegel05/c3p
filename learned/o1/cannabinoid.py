"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by specific structural features such as
    the alkylresorcinol moiety linked to a terpene unit (e.g., in THC and CBD),
    long-chain polyunsaturated fatty acid derivatives linked to ethanolamine
    or glycerol (e.g., in endocannabinoids), and synthetic cannabinoids with
    varied core structures but common pharmacophores.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for phytocannabinoids like THC, CBD (alkylresorcinol linked to a terpene unit)
    phytocannabinoid_pattern = Chem.MolFromSmarts('Oc1cc(ccc1CCCCC)C2=C(C)C=CC=C2C')
    # Include variation in side chain length and saturation
    phytocannabinoid_variations = Chem.MolFromSmarts('Oc1cc(ccc1CCCC)C2=C(C)C=CC=C2C')
    
    # Define pattern for endocannabinoids like anandamide (fatty acid amide linked to ethanolamine)
    endocannabinoid_amide_pattern = Chem.MolFromSmarts('C(=O)NCCO')
    # Define pattern for endocannabinoids like 2-AG (fatty acid ester linked to glycerol)
    endocannabinoid_ester_pattern = Chem.MolFromSmarts('C(=O)OCC(O)CO')

    # Define pattern for synthetic cannabinoids (e.g., CP-55,940)
    synthetic_cannabinoid_pattern = Chem.MolFromSmarts('c1cc(CC(C)(C)C)ccc1C2CCC(CC2)O')

    # Check for phytocannabinoid features
    if mol.HasSubstructMatch(phytocannabinoid_pattern) or mol.HasSubstructMatch(phytocannabinoid_variations):
        return True, "Contains alkylresorcinol moiety linked to a terpene unit characteristic of phytocannabinoids"

    # Check for endocannabinoid features
    if mol.HasSubstructMatch(endocannabinoid_amide_pattern):
        # Verify long-chain fatty acid (at least 16 carbons in chain)
        chain_length = _get_fatty_acid_chain_length(mol, endocannabinoid_amide_pattern)
        if chain_length >= 16:
            return True, "Contains long-chain fatty acid amide linked to ethanolamine characteristic of endocannabinoids"
        else:
            return False, "Amide group found but fatty acid chain too short"

    if mol.HasSubstructMatch(endocannabinoid_ester_pattern):
        # Verify long-chain fatty acid (at least 16 carbons in chain)
        chain_length = _get_fatty_acid_chain_length(mol, endocannabinoid_ester_pattern)
        if chain_length >= 16:
            return True, "Contains long-chain fatty acid ester linked to glycerol characteristic of endocannabinoids"
        else:
            return False, "Ester group found but fatty acid chain too short"

    # Check for synthetic cannabinoid features
    if mol.HasSubstructMatch(synthetic_cannabinoid_pattern):
        return True, "Contains structural features characteristic of synthetic cannabinoids"

    return False, "No characteristic cannabinoid structural features found"

def _get_fatty_acid_chain_length(mol, pattern):
    """
    Helper function to determine the length of the fatty acid chain
    in endocannabinoids.

    Args:
        mol: RDKit Mol object
        pattern: Substructure pattern to match (amide or ester group)

    Returns:
        int: Length of the carbon chain attached to the carbonyl carbon
    """
    chain_lengths = []
    for match in mol.GetSubstructMatches(pattern):
        carbonyl_carbon_idx = match[0]  # Index of the carbonyl carbon
        chain_length = _traverse_chain(mol, carbonyl_carbon_idx)
        chain_lengths.append(chain_length)
    return max(chain_lengths) if chain_lengths else 0

def _traverse_chain(mol, start_idx):
    """
    Traverses a carbon chain starting from the given atom index.

    Args:
        mol: RDKit Mol object
        start_idx: Starting atom index

    Returns:
        int: Length of the carbon chain
    """
    visited = set()
    to_visit = [start_idx]
    chain_length = 0

    while to_visit:
        current_idx = to_visit.pop()
        if current_idx in visited:
            continue
        visited.add(current_idx)
        atom = mol.GetAtomWithIdx(current_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_length += 1
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited]
            to_visit.extend(neighbors)
    return chain_length