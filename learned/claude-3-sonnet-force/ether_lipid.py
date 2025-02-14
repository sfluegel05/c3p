"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:38273 ether lipid

An ether lipid is defined as a lipid similar in structure to a glycerolipid but in which 
one or more of the carbon atoms on glycerol is bonded to an alkyl chain via an ether linkage,
as opposed to the usual ester linkage.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ether oxygen (-O-)
    ether_pattern = Chem.MolFromSmarts("[OX2]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Ether lipids must have at least one ether linkage and no ester linkages
    if not ether_matches or ester_matches:
        return False, "Incorrect linkage pattern for ether lipid"
    
    # Check for long alkyl chains
    alkyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    alkyl_chain_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    if not alkyl_chain_matches:
        return False, "No alkyl chains found"
    
    # Count carbons in the longest alkyl chain
    longest_chain_length = max([len(chain) for chain in alkyl_chain_matches])
    if longest_chain_length < 8:
        return False, "Chains too short for ether lipid"
    
    # Check molecular weight - ether lipids typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for ether lipid"
    
    return True, "Contains at least one ether linkage and long alkyl chains"