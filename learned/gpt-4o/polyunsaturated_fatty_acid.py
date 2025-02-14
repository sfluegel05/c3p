"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group -COOH
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of isolated double bonds in aliphatic chains
    undefined_chain_pattern = Chem.MolFromSmarts("C(=C)C")
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    
    if len(double_bond_matches) <= 1:
        return False, f"Found {len(double_bond_matches)} double bonds, need more than one for polyunsaturation"

    # Confirm that the double bonds are not entirely part of aromatic rings
    if mol.GetNumAromaticRings() > 0:
        return False, "Contains aromatic structures, not typical for fatty acids"
    
    # Check for a sufficiently long aliphatic carbon chain (at least 12)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic())
    if carbon_count < 12:
        return False, f"Chain is too short with only {carbon_count} carbons, typical fatty acids have longer chains"
    
    # Ensure linearity within the chains (no cross-links or complex branching)
    branching_pattern = Chem.MolFromSmarts("CC(C)(C)C")
    if mol.HasSubstructMatch(branching_pattern):
        return False, "Excessive branching interrupts the linear unsaturated chain"
    
    return True, "Contains a carboxylic acid group with more than one isolated aliphatic double bond"