"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a given SMILES string represents a polymer.
    Polymers are characterized by macromolecular chains with repeating substructures,
    typically long carbon chains, and potentially extensive branching.

    Args:
        smiles (str): SMILES string of the chemical entity

    Returns:
        bool: True if the entity is classified as a polymer, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for common repeating units in polymers
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    
    # Identify matches for these patterns
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    amide_matches = len(mol.GetSubstructMatches(amide_pattern))
    ether_matches = len(mol.GetSubstructMatches(ether_pattern))
    
    # Classify as polymer if there are a significant number of repeating units
    repeating_unit_threshold = 3  # Define a threshold for number of repeating units
    if ester_matches >= repeating_unit_threshold or amide_matches >= repeating_unit_threshold or ether_matches >= repeating_unit_threshold:
        return True, "Contains multiple repeating units typical of polymer structure"
    
    # Check for long carbon chains and potential for significant branching
    chain_pattern = Chem.MolFromSmarts("[C]~[C]~[C]~[C]")
    chain_matches = len(mol.GetSubstructMatches(chain_pattern))

    num_c_c_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6)
    num_branching_points = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() > 2)

    # Define molecular weight threshold with caution
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)

    # Enhanced criteria based on chain length and branching
    if (chain_matches >= 3 and num_c_c_bonds > 25 and num_branching_points > 2) or mol_weight > 1000:
        return True, "Has long chains, branching, or high molecular weight, suggesting polymer structure"

    return False, "Does not meet polymer criteria of multiple repeating units or extensive chaining and branching"