"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids are organic aromatic compounds with a structure based on a phenyl ring attached
    to small carbon chains or functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a flexible phenylpropanoid pattern:
    # The main aromatic ring with varying potential linkers (C chains or linking groups)
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
    functional_groups_pattern = Chem.MolFromSmarts("[OX2,CX3,CX4]")
    
    # Check for the presence of at least one phenyl ring
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No aromatic (phenyl) ring found"
        
    # Validate the connection to small aliphatic/functional groups
    if not mol.HasSubstructMatch(functional_groups_pattern):
        return False, "Missing connection to appropriate aliphatic chains or functional groups"

    return True, "Structure consistent with phenylpropanoid, containing phenyl ring and typical linkages"