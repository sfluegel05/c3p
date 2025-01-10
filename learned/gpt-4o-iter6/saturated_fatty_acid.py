"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid should have: 
    - A terminal carboxyl group
    - No carbon-carbon multiple bonds
    - A sufficiently long saturated hydrocarbon chain
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for terminal carboxylic acid group with a long alkyl chain
    carboxylic_pattern = Chem.MolFromSmarts("CC(=O)[O-]")  # Match carboxyl group attached to carbon
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Ensure absence of carbon-carbon multiple bonds
    for bond in mol.GetBonds():
        begin_atom, end_atom = bond.GetBeginAtom(), bond.GetEndAtom()
        if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                return False, "Contains carbon-carbon multiple bond indicating unsaturation"
    
    # Count the number of linearly connected carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbon_atoms)
    if num_carbons < 12:  # Fatty acids generally have longer carbon chains
        return False, "Insufficient carbon length for typical fatty acid"
    
    # Validate the proper hydrocarbon chain pattern with terminal carboxylic acid
    chain_pattern = Chem.MolFromSmarts("CC(=O)[O-]")  # Ensure alkyls form a linear chain ending in carboxyl
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Does not contain sufficiently long saturated carbon chain"

    return True, "Classified as a saturated fatty acid with appropriate structure"