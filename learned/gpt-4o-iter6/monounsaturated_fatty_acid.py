"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    MUFAs have one double or triple bond in the fatty acid chain and  
    a carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Calculate number of carbon-carbon double (or triple) bonds
    double_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    triple_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C#C")))
    
    if (double_bond_count + triple_bond_count) != 1:
        return False, f"Found {double_bond_count} double bonds and {triple_bond_count} triple bonds, need exactly one"

    # Count carbon atoms in the chain (excluding carboxyl group carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() == False)
    if carbon_count < 10:
        return False, "Not enough carbon atoms for a typical fatty acid chain"
    
    return True, "Molecule is a monounsaturated fatty acid (one double or triple bond in the chain with carboxylic group)"