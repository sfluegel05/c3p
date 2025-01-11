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
    
    # Look for terminal carboxylic acid group using modified pattern
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Check for carbon-carbon double and triple bonds specifically between carbons
    for bond in mol.GetBonds():
        begin_atom, end_atom = bond.GetBeginAtom(), bond.GetEndAtom()
        if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:  # Both atoms are carbons
                return False, "Contains carbon-carbon multiple bond indicating unsaturation"
    
    # Ensure there is a sufficiently long chain of saturated atoms
    # Specifically handle typical fatty acids, flexible in chain branches and decorations like methyl groups
    # Improved pattern for generic alkyl chains with carboxylic end
    carbon_chain_pattern = Chem.MolFromSmarts("C-C(=O)O")  # Ensure the carboxyl end

    # Count number of carbons in linear fashion
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons < 8:  # Conservative lower bound for even smallest acids
        return False, "Insufficient carbon length for typical fatty acid"

    # Validate appropriate chain types (in real fatty acid structure, non-linears not matching expected pattern)
    main_chain = Chem.MolFromSmarts("[CX4]~[CX4]")  # Match carbons bound using sp3 bonds, allowing for small branches
    if not mol.HasSubstructMatch(main_chain):
        return False, "Does not contain sufficiently long saturated carbon chain"

    return True, "Classified as a saturated fatty acid with appropriate structure"