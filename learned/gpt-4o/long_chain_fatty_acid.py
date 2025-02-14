"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.

    A long-chain fatty acid is defined as a fatty acid with a carbon chain length from C13 to C22.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Calculate total number of carbons in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Verify the carbon chain length from C13 to C22 (considering COOH group)
    if carbon_count < 13 or carbon_count > 22:
        return False, f"Total carbon count {carbon_count} is not within the long-chain fatty acid range (considering COOH group)"
    
    # Using the calculated exact molecular weight to cross-validate
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if not (200 < mol_weight < 400):  # Typical range for long-chain fatty acids
        return False, "Molecular weight not typical for long-chain fatty acid"
    
    # Check the existence of a large linear sequence of carbons
    large_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(large_chain_pattern):
        return False, "Insufficient chain length or incorrect structure for fatty acid"
    
    return True, f"Valid long-chain fatty acid with {carbon_count} carbon atoms"