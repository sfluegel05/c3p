"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: CHEBI:82198 wax
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdchem

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is defined as an organic compound or mixture of compounds
    that is composed of long-chain molecules and is malleable at ambient temperatures.
    Commonly, waxes are esters formed from long-chain fatty acids and long-chain alcohols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Attempt to parse SMILES with error handling
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Chem.rdchem.KekulizeException as e:
        return False, f"SMILES parsing error: {str(e)}"
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify ester groups using SMARTS pattern
    ester_smarts = '[CX3](=O)[OX2H0][#6]'
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Define threshold for long-chain (number of carbons)
    chain_length_threshold = 12

    # For each ester group, split the molecule and analyze fragments
    for match in ester_matches:
        carbonyl_c_idx = match[0]
        ester_o_idx = match[2]

        # Break the ester bond (between carbonyl carbon and oxygen)
        rw_mol = Chem.RWMol(mol)
        rw_mol.BeginBatchEdit()
        rw_mol.RemoveBond(carbonyl_c_idx, ester_o_idx)
        rw_mol.CommitBatchEdit()

        # Get the fragments resulting from bond breakage
        frags = Chem.GetMolFrags(rw_mol, asMols=True, sanitizeFrags=True)

        # Initialize counts for long chains
        long_chain_count = 0

        # Analyze each fragment
        for frag in frags:
            # Calculate the longest aliphatic carbon chain
            longest_chain = get_longest_aliphatic_chain(frag)
            
            if longest_chain >= chain_length_threshold:
                long_chain_count += 1

        # If both fragments have long aliphatic chains, it's a wax
        if long_chain_count >= 2:
            return True, f"Ester with two long aliphatic chains found (long chains: {long_chain_count})"

    return False, "No ester with two long aliphatic chains found"

def get_longest_aliphatic_chain(mol):
    """
    Calculates the length of the longest continuous aliphatic carbon chain in a molecule.

    Args:
        mol (rdkit.Chem.Mol): Molecule to analyze

    Returns:
        int: Length of the longest aliphatic carbon chain
    """
    # Find all carbon atoms
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]

    # Build a graph of carbon-carbon connections
    paths = []
    for carbon_idx in carbon_idxs:
        for other_carbon_idx in carbon_idxs:
            if carbon_idx < other_carbon_idx:
                # Find all paths between two carbon atoms
                all_paths = Chem.rdmolops.GetAllPaths(mol, carbon_idx, other_carbon_idx)
                for path in all_paths:
                    # Check if all atoms in path are carbons
                    if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in path):
                        # Check if bonds are single or double (allow for unsaturation)
                        if all(mol.GetBondBetweenAtoms(path[i], path[i+1]).GetBondType() in [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE] for i in range(len(path)-1)):
                            paths.append(path)

    # Find the longest path
    if paths:
        longest_path_length = max(len(path) for path in paths)
    else:
        longest_path_length = 0

    return longest_path_length