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

    # Count the number of double bonds specifically in aliphatic chains
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) <= 1:
        return False, f"Found {len(double_bond_matches)} double bonds, need more than one for polyunsaturation"

    # Confirm that the double bonds are not part of aromatic structures
    if mol.GetNumAromaticRings() > 0:
        return False, "Contains aromatic structures, not typical for fatty acids"
    
    # Check for a sufficiently long aliphatic carbon chain (at least 12)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() == False)
    if carbon_count < 12:
        return False, f"Chain is too short with only {carbon_count} carbons, typical fatty acids have longer chains"

    # Check if double bonds belong to one chain/sequential without being terminated by branches
    saturated_chain = Chem.MolFromSmarts("CC(C)=C")  # Example of branching interrupt
    if mol.HasSubstructMatch(saturated_chain):
        return False, "Branching interrupts the linear unsaturated chain"
    
    return True, "Contains a carboxylic acid group with more than one double bond in a non-aromatic, linear chain"