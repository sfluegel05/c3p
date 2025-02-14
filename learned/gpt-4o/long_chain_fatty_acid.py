"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

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

    # Calculate the molecular weight
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    # Adjust molecular weight range for validation
    if not (180 < mol_weight < 450):
        return False, "Molecular weight not typical for long-chain fatty acid"
    
    # Calculate total number of carbons in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Check for total number of carbon chains compatible with a C13-C22 length, potentially excluding certain complex functional groups
    # Consider flexible identification of the carbon chain for possible alkenes, alkynes, and partial branching without ring systems
    large_chain_pattern = Chem.MolFromSmarts("C(-C)~C~C~C~C~C~C~C~C~C")
    if not mol.HasSubstructMatch(large_chain_pattern):
        return False, "Insufficient chain length or incorrect fatty acid structure detected"
    
    # Reassess the carbon count for the core chain
    if carbon_count < 13 or carbon_count > 22:
        return False, f"Total carbon count {carbon_count} is not within the long-chain fatty acid range"

    return True, f"Valid long-chain fatty acid with {carbon_count} carbon atoms"