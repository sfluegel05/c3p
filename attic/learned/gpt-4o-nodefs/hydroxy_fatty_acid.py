"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid has at least one hydroxyl group (-OH) and a carboxylic acid group (-COOH) with a hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for at least one hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"
    
    # Check hydrocarbon chain length (at least 8 carbon atoms)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 8:
        return False, f"Insufficient carbon atoms for a fatty acid chain, found {carbon_count}"
    
    # Additional checks for hydroxy fatty acids characteristics can be added as needed
    
    return True, "Contains the required hydroxyl and carboxylic acid groups with sufficient hydrocarbon chain"

# Example usage to test the function
example_smiles = "CC(CO)CCCC(C)CC(O)=O"  # Example: omega-hydroxyphytanic acid
result, reason = is_hydroxy_fatty_acid(example_smiles)
print(f"Is hydroxy fatty acid: {result}, Reason: {reason}")