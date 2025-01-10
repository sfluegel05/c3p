"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: Clavulone (a class of esterified prostanoids obtained from marine corals)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    Clavulones are esterified prostanoids with a cyclopentane ring, hydroxyl group, and ester groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cyclopentane ring with a ketone group (C=O) and hydroxyl group (OH)
    # Use a more flexible pattern to account for substitutions and stereochemistry
    clavulone_core_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1(=O)[OH]")
    if not mol.HasSubstructMatch(clavulone_core_pattern):
        return False, "No cyclopentane ring with ketone and hydroxyl group found"

    # Check for ester groups (acetates or similar) with specific patterns
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found"

    # Check for long hydrocarbon chains with double bonds (C=C) specific to clavulones
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    hydrocarbon_chain_matches = mol.GetSubstructMatches(hydrocarbon_chain_pattern)
    if len(hydrocarbon_chain_matches) < 1:
        return False, "No long hydrocarbon chains with double bonds found"

    # Check for the presence of at least one double bond in the hydrocarbon chain
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 1:
        return False, "No double bonds found in hydrocarbon chains"

    # Check for specific arrangement of functional groups around the cyclopentane ring
    functional_group_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1(=O)[OH][CX3](=[OX1])[CX4]")
    if not mol.HasSubstructMatch(functional_group_pattern):
        return False, "No specific arrangement of functional groups around the cyclopentane ring"

    return True, "Contains cyclopentane ring with ketone, hydroxyl, ester groups, and long hydrocarbon chains with double bonds"