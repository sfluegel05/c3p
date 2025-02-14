"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: Epoxy Fatty Acid
"""
from rdkit import Chem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a heterocyclic fatty acid containing an epoxide ring as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of aromatic rings (should not have any)
    if mol.GetNumAromaticRings() > 0:
        return False, "Molecule contains aromatic rings, which is not characteristic of fatty acids"

    # Check for carboxylic acid group (-C(=O)OH) at the end of the chain
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    ca_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not ca_matches:
        return False, "No carboxylic acid group found"

    # Identify the longest aliphatic chain ending with the carboxylic acid
    # Define a fatty acid chain as a chain of carbons without branching (degree <=2)
    fatty_acid_chain = None
    max_chain_length = 0
    for ca_match in ca_matches:
        ca_atom_idx = ca_match[0]  # Carbon atom of the carboxylic group
        paths = Chem.rdmolops.FindAllPathsOfLengthN(mol, 2, useBonds=False)
        # Use BFS to find the longest chain ending at the carboxylic carbon
        visited = set()
        queue = [(ca_atom_idx, [ca_atom_idx])]
        while queue:
            current_idx, path = queue.pop(0)
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            if atom.GetDegree() > 2 and current_idx != ca_atom_idx:
                continue  # Skip branching
            neighbors = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 6]
            for neighbor in neighbors:
                if neighbor not in visited:
                    new_path = path + [neighbor]
                    queue.append((neighbor, new_path))
                    if len(new_path) > max_chain_length:
                        max_chain_length = len(new_path)
                        fatty_acid_chain = new_path

    if fatty_acid_chain is None or max_chain_length < 12:
        return False, f"Molecule does not have a long enough aliphatic chain (found length {max_chain_length})"

    # Check for epoxide ring (three-membered ring with one oxygen and two carbons)
    epoxide_pattern = Chem.MolFromSmarts("[C;R1]1[C;R1][O;R1]1")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if not epoxide_matches:
        return False, "No epoxide ring found"

    # Ensure that the epoxide ring is part of the aliphatic chain
    epoxide_in_chain = False
    chain_atom_indices = set(fatty_acid_chain)
    for epoxide_match in epoxide_matches:
        if any(idx in chain_atom_indices for idx in epoxide_match):
            epoxide_in_chain = True
            break

    if not epoxide_in_chain:
        return False, "Epoxide ring is not part of the aliphatic chain"

    return True, "Molecule is an epoxy fatty acid: contains a long aliphatic chain with a terminal carboxylic acid and an epoxide ring"