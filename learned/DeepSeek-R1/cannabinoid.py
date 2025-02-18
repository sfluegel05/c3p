"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: CHEBI:64211 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by oxygen-containing functional groups or heterocyclic rings,
    and include endocannabinoids like anandamide or phytocannabinoids like THC.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for ethanolamide group (N-linked to CO and ethylene glycol)
    ethanolamide_pattern = Chem.MolFromSmarts("[NH]CCO")
    if mol.HasSubstructMatch(ethanolamide_pattern):
        return True, "Contains ethanolamide group characteristic of endocannabinoids"
    
    # Check for glycerol ester structure (like 2-AG)
    glycerol_ester_pattern = Chem.MolFromSmarts("[O;R][C](=O)CCCCCC")
    if mol.HasSubstructMatch(glycerol_ester_pattern):
        return True, "Contains glycerol ester with long-chain fatty acid"
    
    # Check for dibenzopyran ring system (characteristic of THC-like structures)
    thc_core = Chem.MolFromSmarts("C12CCC(C)([C@@H]3C=CC(OC1)=CC3)CC2")
    if mol.HasSubstructMatch(thc_core):
        return True, "Contains dibenzopyran cannabinoid core structure"
    
    # Check for oxygen in heterocyclic rings (e.g., tetrahydrofuran, epoxide)
    oxygen_in_ring = any(atom.GetAtomicNum() == 16 and atom.IsInRing() for atom in mol.GetAtoms())
    if oxygen_in_ring:
        return True, "Contains oxygen in a heterocyclic ring"
    
    # General check for oxygen-containing functional groups (ester, ether, etc.)
    oxygen_fgs = Chem.MolFromSmarts("[O;!H0]")
    if mol.HasSubstructMatch(oxygen_fgs):
        return True, "Contains oxygen functional group typical of cannabinoids"
    
    return False, "Does not match known cannabinoid structural patterns"