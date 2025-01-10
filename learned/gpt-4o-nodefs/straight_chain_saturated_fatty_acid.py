"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the carboxylic acid group (-C(=O)O)
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylic acid group found"
    
    # Identify the largest carbon chain
    carbon_chains = mol.GetSubstructMatches(Chem.MolFromSmarts("C"))
    if not carbon_chains:
        return False, "No carbon chain found"

    # Check for saturation (no double or triple bonds in the carbon chain)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            if any([bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE) for bond in atom.GetBonds()]):
                return False, "Carbon chain contains double or triple bonds"

    # Check for branching (more than two single bonds for any carbon)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            n_single_bonds = sum(1 for bond in atom.GetBonds() if bond.GetBondType() == Chem.BondType.SINGLE)
            if n_single_bonds > 2:
                return False, "Carbon chain is branched"
    
    return True, "Molecule is a straight-chain saturated fatty acid"

# Examples
example_smiles = "CCCCCCC(O)=O"  # Heptanoic acid
result, reason = is_straight_chain_saturated_fatty_acid(example_smiles)
print(f"{example_smiles}: {result} ({reason})")