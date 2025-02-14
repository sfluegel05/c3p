"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: trienoic fatty acid
Definition: Any polyunsaturated fatty acid that contains three double bonds.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a fatty acid with a terminal carboxylic acid group and at least three
    carbon-carbon double bonds in its longest unbranched aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Neutralize charges
    Chem.SanitizeMol(mol)

    # Look for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    carbox_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    if not carbox_matches:
        return False, "No terminal carboxylic acid group found"

    # Identify the carboxylic carbon atom
    carbox_carbons = [match[0] for match in carbox_matches]

    # Find all aliphatic carbon chains ending with carboxylic acid group
    chains = []
    for carbox_c in carbox_carbons:
        paths = Chem.rdmolops.FindAllPathsOfLengthN(mol, mol.GetNumAtoms(), useBonds=False)
        # Filter paths that start from the carboxylic carbon
        for path in paths:
            if path[-1] == carbox_c:
                chains.append(path)

    # If no chains found, consider single bonds from carboxylic carbon
    if not chains:
        chains = Chem.rdmolops.FindAllPathsOfLengthN(mol, mol.GetNumAtoms(), useBonds=False, useHs=True)

    # Filter for the longest chain
    longest_chain = []
    for path in chains:
        # Check if path consists only of carbons (except for the last carboxylic oxygen atoms)
        atoms_in_path = [mol.GetAtomWithIdx(idx) for idx in path]
        if all(atom.GetAtomicNum() == 6 for atom in atoms_in_path[:-1]):
            if len(path) > len(longest_chain):
                longest_chain = path

    if not longest_chain:
        return False, "No suitable aliphatic chain found"

    # Extract the substructure corresponding to the longest chain
    chain_atoms = set(longest_chain)
    chain_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in chain_atoms and bond.GetEndAtomIdx() in chain_atoms:
            chain_bonds.append(bond)

    # Count the number of carbon-carbon double bonds in the chain
    num_double_bonds = 0
    for bond in chain_bonds:
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if (bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and
            begin_atom.GetAtomicNum() == 6 and
            end_atom.GetAtomicNum() == 6):
            num_double_bonds += 1

    if num_double_bonds < 3:
        return False, f"Contains {num_double_bonds} carbon-carbon double bonds in the chain, needs at least 3"

    # Check for branching (degree > 2 for any carbon except the carboxylic carbon)
    for idx in longest_chain[:-1]:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetDegree() > 2:
            return False, "Chain is branched, fatty acids are typically unbranched"

    # Count carbons in the chain
    c_count = sum(1 for idx in longest_chain if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if c_count < 12:
        return False, f"Chain contains {c_count} carbon atoms, which is too short for a typical fatty acid"

    return True, "Molecule is a trienoic fatty acid with a terminal carboxylic acid group and at least three C=C double bonds in the aliphatic chain"