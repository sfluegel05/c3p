"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: CHEBI:37577 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed from the condensation of the thiol group of coenzyme A
    with the carboxy group of any carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a simplified Coenzyme A SMARTS pattern
    # This pattern captures key features of CoA: adenine ring, ribose, phosphates, pantetheine chain ending with thiol
    coa_smarts = '[#7]-[#6](=O)-[#6]-[#6]-[#7]-[#6](=O)-[#6]([#8])-[#6]([#6])([#6])-[#6]-[#8]-[#6]-[#8]-[#8]-[#6]-[#8]-c1nc2c(n1)nc(nc2N)[#7]-[#6]-1-[#6]([#8])-[#6]([#8])-[#6]-[#8]-1'

    coa_mol = Chem.MolFromSmarts(coa_smarts)
    if coa_mol is None:
        return False, "Failed to parse Coenzyme A SMARTS pattern"
    
    # Check for Coenzyme A moiety in the molecule, ignoring chirality
    if not mol.HasSubstructMatch(coa_mol, useChirality=False):
        return False, "Coenzyme A moiety not found"
    
    # Define SMARTS pattern for the thioester linkage connected to the CoA sulfur
    thioester_smarts = '[#6][C](=O)[S]-[#6]'
    thioester_mol = Chem.MolFromSmarts(thioester_smarts)
    if thioester_mol is None:
        return False, "Failed to parse thioester SMARTS pattern"
    
    # Find matches for thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_mol, useChirality=False)
    if not thioester_matches:
        return False, "Thioester linkage not found"
    
    # Ensure the sulfur in the thioester is connected to the CoA moiety
    for match in thioester_matches:
        sulfur_idx = match[2]  # Index of the sulfur atom in the match
        # Check if sulfur atom is part of the CoA moiety
        atom_map = {}
        if mol.HasSubstructMatch(coa_mol, uniChemistry.py.view_map(0s)), atomMap=atom_map
                                    , useChirality=False):
            if sulfur_idx in atom_map.values():
                return True, "Contains Coenzyme A moiety linked via thioester bond to an acyl group"
    
    return False, "Thioester linkage not connected to CoA sulfur atom"