"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    MUFAs have one double or triple bond in the longest carbon chain, typically with a carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find the longest carbon chain
    longest_chain = Chem.rdmolops.GetLongestChain(mol)
    if len(longest_chain) == 0:
        return False, "No valid carbon chain found"

    # Extract substructure matching the longest chain
    chain_atoms = [mol.GetAtomWithIdx(idx) for idx in longest_chain]
    chain_bonds = [mol.GetBondBetweenAtoms(chain_atoms[i].GetIdx(), chain_atoms[i + 1].GetIdx())
                   for i in range(len(chain_atoms) - 1)]
    
    # Count the number of double or triple bonds in the longest chain
    doub_triple_bond_count = sum(bond.GetBondTypeAsDouble() in [2.0, 3.0] for bond in chain_bonds)
    if doub_triple_bond_count != 1:
        return False, f"Found {doub_triple_bond_count} double/triple bonds in the longest chain, require exactly one"

    # Verify that all atoms in the longest chain are carbon (except functional groups)
    non_carbon_atoms = [atom for atom in chain_atoms if atom.GetAtomicNum() != 6 and atom.GetIdx() not in longest_chain]
    if non_carbon_atoms:
        return False, "Contains non-carbon atoms in the main chain, other than functional groups"

    # Check for the presence of a carboxylic acid group at the terminal
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"

    return True, "Molecule is a monounsaturated fatty acid with one double/triple bond in the main chain"