"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: long-chain fatty alcohol (C13-C22)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol (C13-C22 chain) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Function to compute the longest carbon chain length
    def get_longest_carbon_chain(mol):
        # Create a copy with only carbons and single bonds
        emol = Chem.RWMol(mol)
        # Remove non-carbon atoms
        atoms_to_remove = [atom.GetIdx() for atom in emol.GetAtoms() if atom.GetAtomicNum() != 6]
        for idx in sorted(atoms_to_remove, reverse=True):
            emol.RemoveAtom(idx)
        # Convert all bonds to single to consider only single bonds for chain length
        for bond in emol.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
        mol_carbons = emol.GetMol()
        return rdMolDescriptors.CalcLongestChainLength(mol_carbons)

    longest_chain = get_longest_carbon_chain(mol)

    # Check chain length is between 13 and 22
    if longest_chain < 13 or longest_chain > 22:
        return False, f"Longest carbon chain is {longest_chain}, not in 13-22"

    return True, f"Long-chain fatty alcohol with C{longest_chain} chain"