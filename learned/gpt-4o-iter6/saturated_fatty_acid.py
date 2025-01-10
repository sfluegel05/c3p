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
    
    # Look for carboxylic acid group (COOH)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Check for carbon-carbon double and triple bonds
    if any(bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE] for bond in mol.GetBonds()):
        return False, "Contains carbon-carbon multiple bond indicating unsaturation"
    
    # Ensure there is a sufficiently long chain of sp3 carbons
    # Here, target at least 6 carbons in a row, allowing for typical saturated chains
    carbon_chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]") 
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Does not contain sufficiently long saturated carbon chain"

    # Account for simple methyl or hydroxyl branches
    # Ensure no complex or non-linear backbone structures that confuse classification
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in {1, 6, 8}:  # Elements other than C, O, H are uncommon in simple acids
            return False, "Contains unusual side groups not typical of simple saturated fatty acids"

    return True, "Classified as a saturated fatty acid with appropriate structure"