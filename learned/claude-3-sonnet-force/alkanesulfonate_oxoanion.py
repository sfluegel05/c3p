"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: CHEBI:52526 alkanesulfonate oxoanion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is defined as having a sulfonate group (-SO3-) attached to a carbon chain,
    where the carbon at position 1 can be attached to various R groups (hydrogens, carbon chains, or other groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sulfonate group pattern
    sulfonate_pattern = Chem.MolFromSmarts("[S+3(=O)(-[O-])(-[O-])]")
    sulfonate_matches = mol.GetSubstructMatches(sulfonate_pattern)
    
    if len(sulfonate_matches) == 0:
        return False, "No sulfonate group found"
    
    # Look for carbon chain attached to the sulfonate group
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3][CX4,CX3]*[S+3(=O)(-[O-])(-[O-])]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern, sulfonate_matches)
    
    if len(carbon_chain_matches) > 0:
        return True, "Molecule contains an alkanesulfonate oxoanion group"
    else:
        return False, "No alkanesulfonate oxoanion group found (missing carbon chain)"