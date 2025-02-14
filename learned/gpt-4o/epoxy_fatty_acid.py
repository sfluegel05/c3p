"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid contains an epoxide ring as part of a long aliphatic chain ending in a carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Confirm presence of an epoxide group within an aliphatic chain (avoid terminal position)
    epoxide_pattern = Chem.MolFromSmarts("C[C@H]1O[C@@H]1C")  # More specific epoxide, part of a chain
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found or not part of an aliphatic chain"
    
    # Check for the presence of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Verify the chain length - minimum 12 carbon atoms
    aliphatic_chain_pattern = Chem.MolFromSmarts("C[CH2]N*"*6)  # Minimum 12 carbons counted
    if not mol.HasSubstructMatch(aliphatic_chain_pattern):
        return False, "Aliphatic chain inadequate to be a fatty acid"
    
    # Ensure exclusion of large ring systems, ensure it is primarily a chain
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:  # The only ring should be the epoxide
        return False, "Contains large ring systems uncharacteristic of fatty acids"

    return True, "Contains epoxide ring as part of a long aliphatic chain ending in a carboxyl group"