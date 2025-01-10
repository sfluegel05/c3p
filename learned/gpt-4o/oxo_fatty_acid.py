"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid contains at least one aldehydic or ketonic group in addition to a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a carboxylic acid group -C(=O)O
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for aldehydic group -[CX3H1]=O
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1]=O")
    
    # Look for ketonic group -[CX3](=O)[#6]
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    
    # Determine presence of aldehydic or ketonic groups
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    has_ketone = mol.HasSubstructMatch(ketone_pattern)

    # If neither is present, it is not an oxo fatty acid
    if not has_aldehyde and not has_ketone:
        return False, "No aldehydic or ketonic group found"
    
    # Ensuring the molecule is predominantly a linear chain (characteristic of fatty acids)
    num_carbon_chains = 0
    chain_pattern = Chem.MolFromSmarts("CCCC")
    if mol.HasSubstructMatch(chain_pattern):
        num_carbon_chains = len(mol.GetSubstructMatches(chain_pattern))

    # Check total number of carbon atoms
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())

    # Typical fatty acids have a linear backbone, count should be reasonable
    if num_carbon_chains < 2 or num_carbons < 5:
        return False, "Molecule does not match typical oxo fatty acid structure"

    return True, "Contains both carboxylic acid group and aldehydic or ketonic group, characteristic of an oxo fatty acid"