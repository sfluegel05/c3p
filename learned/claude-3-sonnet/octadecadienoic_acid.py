"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid
Definition: Any straight-chain, C18 polyunsaturated fatty acid having two C=C double bonds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbons in main chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 18:
        return False, f"Contains {carbon_count} carbons, must be exactly 18"
    
    # Count double bonds (exclude the one in carboxylic acid)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    
    # Need exactly 2 C=C double bonds (not counting the C=O in carboxylic acid)
    if double_bonds != 2:
        return False, f"Contains {double_bonds} C=C double bonds, must be exactly 2"
    
    # Check if it's a straight chain (main carbon backbone)
    # First, get the carboxyl carbon
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "Could not identify carboxyl group"
    
    # Count carbons in longest chain from carboxyl group
    carboxyl_carbon = carboxyl_matches[0][0]  # First atom of first match
    longest_chain = Chem.rdmolops.GetShortestPath(mol, carboxyl_carbon, -1)
    max_chain_length = max(len(path) for path in longest_chain if path)
    
    # The longest chain should include most of the carbons
    # (allowing for some substitutions)
    if max_chain_length < 15:  # Using 15 as threshold to allow for some branching
        return False, "Not a straight chain structure"
    
    # Check for cyclic structures
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Contains rings, must be acyclic"
        
    return True, "C18 straight-chain fatty acid with 2 C=C double bonds"