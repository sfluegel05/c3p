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
    
    # Linear chain pattern of aliphatic backbone: greater than 4 carbon chain
    long_chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No sufficiently long carbon chain typical of fatty acids found"
    
    # Verify types of bonded atoms around the carbonyl groups
    # Ensuring at least one terminal carbonyl group indicative of long chain
    terminal_ketone_pattern = Chem.MolFromSmarts("[CH2]!@C(=O)[CH2]")
    if not mol.HasSubstructMatch(terminal_ketone_pattern):
        return False, "No terminal ketone group along carbon chain"

    # This condition will also help to eliminate some macromolecules or polysaccharides incorrectly classified
    
    return True, "Contains both carboxylic acid group and aldehydic or ketonic group, characteristic of an oxo fatty acid"