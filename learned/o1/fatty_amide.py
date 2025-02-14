"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI:24038 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define amide functional group pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    
    # Find all amide groups in the molecule
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amides = len(amide_matches)
    
    if num_amides == 0:
        return False, "No amide functional group found"
    elif num_amides > 1:
        return False, f"Found {num_amides} amide groups, expecting only one for monocarboxylic acid amide"
    
    # Get the carbonyl carbon atom index of the amide
    amide_match = amide_matches[0]
    carbonyl_c_idx = amide_match[0]
    carbonyl_n_idx = amide_match[2]
    
    # Traverse the chain attached to the carbonyl carbon (excluding the carbonyl oxygen and nitrogen)
    chain_atoms = set()
    atoms_to_visit = [carbonyl_c_idx]
    visited_atoms = set()
    
    while atoms_to_visit:
        atom_idx = atoms_to_visit.pop()
        if atom_idx in visited_atoms:
            continue
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            continue  # Skip if not carbon
        if atom_idx == carbonyl_c_idx:
            # Skip the carbonyl oxygen and nitrogen
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 
                         and nbr.GetIdx() != amide_match[1] and nbr.GetIdx() != carbonyl_n_idx]
        else:
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 
                         and nbr.GetIdx() not in visited_atoms]
        chain_atoms.add(atom_idx)
        atoms_to_visit.extend(neighbors)
    
    chain_length = len(chain_atoms)
    
    if chain_length < 6:
        return False, f"Chain length is {chain_length}, which is too short for a fatty acid chain"
    
    # Check if the chain is unbranched (no side chains)
    is_branched = False
    for atom_idx in chain_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        carbon_neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 2:
            is_branched = True
            break
    if is_branched:
        return False, "The fatty acid chain is branched, fatty acids are typically unbranched"
    
    # Check for presence of functional groups in the chain
    has_functional_groups = False
    for atom_idx in chain_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            has_functional_groups = True
            break
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() not in [1,6]:
                has_functional_groups = True
                break
    if has_functional_groups:
        return False, "Functional groups detected in the fatty acid chain"
    
    return True, "Molecule is a fatty amide with a long unbranched hydrocarbon chain"

# Example usage:
# smiles = "CCCCCCCCCCCCCCCC(=O)NCCO"  # N-palmitoyl ethanolamide
# result, reason = is_fatty_amide(smiles)
# print(result, reason)