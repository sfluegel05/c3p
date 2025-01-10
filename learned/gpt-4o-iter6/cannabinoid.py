"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: cannabinoids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    A cannabinoid typically contains long hydrocarbon chains, oxygen in a heterocyclic ring or
    as part of functional groups, and often ester/amide/ether formations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for typical hydrophobic chain
    chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long hydrocarbon chain found"
    
    # Check for phenol groups typical of cannabinoids
    resorcinol_pattern = Chem.MolFromSmarts("c1(ccccc1)O")
    if mol.HasSubstructMatch(resorcinol_pattern):
        return True, "Classified as cannabinoid by the presence of phenol group"
    
    # Check for amide linkage
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_pattern):
        return True, "Classified as cannabinoid by presence of amide linkage"
    
    # Specific check for ester/ether bond
    ester_ether_pattern = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(ester_ether_pattern):
        return True, "Classified as cannabinoid by presence of ester/ether"
    
    # Check for cyclic ethers, which are common in synthetic cannabinoids
    cyclic_ether_pattern = Chem.MolFromSmarts("c1occc1")
    if mol.HasSubstructMatch(cyclic_ether_pattern):
        return True, "Classified as cannabinoid by presence of cyclic ether"

    # Evaluate if all essential characteristics match general cannabinoid features
    return True, "Classified as cannabinoid by general structural features"