"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: CHEBI:38867 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    A bisbenzylisoquinoline alkaloid consists of two benzylisoquinoline units linked by ether bridges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core structure of a benzylisoquinoline unit
    benzylisoquinoline_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CH2]-[c]2[c][c][c][c][c]2-[N]")
    if not mol.HasSubstructMatch(benzylisoquinoline_pattern):
        return False, "No benzylisoquinoline unit found"

    # Look for two benzylisoquinoline units
    matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(matches) < 2:
        return False, f"Found {len(matches)} benzylisoquinoline units, need at least 2"

    # Look for ether bridges (O-C) between the units
    ether_bridge_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    ether_bridge_matches = mol.GetSubstructMatches(ether_bridge_pattern)
    if len(ether_bridge_matches) < 1:
        return False, "No ether bridges found between benzylisoquinoline units"

    # Check for additional bridging (carbon-carbon or methylenedioxy)
    carbon_carbon_bridge_pattern = Chem.MolFromSmarts("[CX4][CX4]")
    methylenedioxy_pattern = Chem.MolFromSmarts("[OX2][CX4][OX2]")
    
    has_carbon_carbon_bridge = mol.HasSubstructMatch(carbon_carbon_bridge_pattern)
    has_methylenedioxy_bridge = mol.HasSubstructMatch(methylenedioxy_pattern)
    
    if not (has_carbon_carbon_bridge or has_methylenedioxy_bridge):
        return False, "No additional bridging (carbon-carbon or methylenedioxy) found"

    # Check molecular weight - bisbenzylisoquinoline alkaloids typically have high molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for bisbenzylisoquinoline alkaloid"

    # Count nitrogen atoms - should have at least 2 (one per benzylisoquinoline unit)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, "Too few nitrogen atoms for bisbenzylisoquinoline alkaloid"

    return True, "Contains two benzylisoquinoline units linked by ether bridges with additional bridging"