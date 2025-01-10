"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    These compounds are derived from tryptophan and typically secologanin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic nitrogen (required for alkaloid)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count == 0:
        return False, "No nitrogen atoms found - required for alkaloid"

    # Look for indole or modified indole core
    indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")  # Basic indole
    modified_indole_pattern = Chem.MolFromSmarts("c1ccc2nccc2c1")  # Modified indole
    if not (mol.HasSubstructMatch(indole_pattern) or mol.HasSubstructMatch(modified_indole_pattern)):
        return False, "No indole or modified indole core found"

    # Check for complex polycyclic structure typical of MIAs
    ring_info = rdMolDescriptors.CalcMolFormula(mol)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 4:
        return False, f"Too few rings ({ring_count}) for MIA structure"

    # Check carbon count (typically >18 for MIA skeleton)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:
        return False, f"Too few carbons ({c_count}) for MIA structure"

    # Check for presence of bridged bicyclic systems common in MIAs
    bridged_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C][C][C]1[C]2")
    if not mol.HasSubstructMatch(bridged_pattern):
        return False, "Missing typical bridged ring systems"

    # Look for typical MIA features:
    # - Often contains methyl ester
    # - Often contains ethyl group
    # - Often contains additional oxygenation
    methyl_ester = Chem.MolFromSmarts("CC(=O)O")
    ethyl = Chem.MolFromSmarts("CC")
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    features = []
    if mol.HasSubstructMatch(methyl_ester):
        features.append("methyl ester")
    if mol.HasSubstructMatch(ethyl):
        features.append("ethyl group")
    if oxygen_count > 0:
        features.append(f"{oxygen_count} oxygen atoms")

    # Check for sp2 carbons (common in terpene portion)
    sp2_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[C]=[C,N,O]")))
    if sp2_carbons == 0:
        return False, "No sp2 carbons found - unusual for MIA"

    # Calculate complexity score
    complexity_contributors = [
        ring_count * 2,  # Rings contribute to complexity
        n_count * 3,     # Nitrogens contribute to complexity
        oxygen_count * 2, # Oxygens contribute to complexity
        sp2_carbons,     # sp2 carbons contribute to complexity
    ]
    complexity_score = sum(complexity_contributors)
    
    if complexity_score < 15:
        return False, f"Complexity score ({complexity_score}) too low for typical MIA"

    feature_str = ", ".join(features) if features else "typical structural features"
    return True, f"Contains indole core, complex polycyclic system ({ring_count} rings), and {feature_str}"