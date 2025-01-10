"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a given SMILES string represents a polymer via enhanced structural 
    complexity detection, diverse repeating unit patterns, and molecular weight considerations.

    Args:
        smiles (str): SMILES string of the chemical entity

    Returns:
        bool: True if the entity is classified as a polymer, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Define additional SMARTS patterns that might indicate polymer-like structures
    aromatic_chain_pattern = Chem.MolFromSmarts("c-c")
    conjugated_system_pattern = Chem.MolFromSmarts("C=C-C=C")
    
    # Identify matches for these patterns
    aromatic_chain_matches = len(mol.GetSubstructMatches(aromatic_chain_pattern))
    conjugated_system_matches = len(mol.GetSubstructMatches(conjugated_system_pattern))
    
    # Count functional groups that can form repeat units
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    amide_matches = len(mol.GetSubstructMatches(amide_pattern))
    ether_matches = len(mol.GetSubstructMatches(ether_pattern))
    
    # Use RDKit descriptors for complexity and size
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    total_rings = Chem.rdMolDescriptors.CalcNumRings(mol)
    
    # Adapt threshold based on mol_weight and structural matches
    repeating_unit_threshold = 3  # This can be adjusted
    
    if (ester_matches >= repeating_unit_threshold or amide_matches >= repeating_unit_threshold or
        ether_matches >= repeating_unit_threshold or aromatic_chain_matches > 5 or
        conjugated_system_matches > 3):
        if (mol_weight > 1000 or rotatable_bonds > 20 or total_rings > 10):
            return (True, "Exhibits multiple repeating units or extensive molecular complexity typical of polymers")
    
    return (False, "Does not meet polymer criteria of multiple repeating units or extensive molecular complexity")