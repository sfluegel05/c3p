"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:136568 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is any fatty acid with a hydroxy functional group in the alpha- or 2-position.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for 2-hydroxy fatty acid
    # Carboxylic acid group connected to an alpha carbon with hydroxy group,
    # followed by an aliphatic chain (variable length)
    pattern = Chem.MolFromSmarts("""
    [C;X3](=O)[O;H1,H0-]   # Carboxylic acid group
    [C;H1]                 # Alpha carbon with one hydrogen
    (
        [O;H1]             # Attached hydroxy group
    )
    [C;H2][C;H2]          # At least two CH2 groups (minimum chain length)
    """)
    
    if not pattern:
        return False, "Invalid SMARTS pattern"

    # Search for matches in the molecule
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No matching 2-hydroxy fatty acid substructure found"

    # Verify each match
    for match in matches:
        carboxy_c_idx = match[0]
        alpha_c_idx = match[2]  # Index of alpha carbon
        hydroxy_o_idx = match[3]

        alpha_c = mol.GetAtomWithIdx(alpha_c_idx)

        # Ensure alpha carbon is only connected to carboxylic carbon, hydroxy oxygen, and next carbon
        neighbors = alpha_c.GetNeighbors()
        if len(neighbors) != 3:
            continue  # More than expected neighbors, not a simple fatty acid

        # Check that the chain is aliphatic and sufficiently long
        chain_length = 0
        current_atom = neighbors[-1]  # Get the carbon after the alpha carbon
        visited_atoms = set([alpha_c_idx])

        while current_atom.GetAtomicNum() == 6 and current_atom.GetDegree() == 2:
            # Carbon atoms in the chain should be sp3 hybridized (single bonds)
            bonds = [bond for bond in current_atom.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.SINGLE]
            if len(bonds) != 2:
                break  # Not a simple chain
            chain_length += 1
            visited_atoms.add(current_atom.GetIdx())

            # Move to the next carbon in the chain
            next_atoms = [nbr for nbr in current_atom.GetNeighbors()
                          if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited_atoms]
            if not next_atoms:
                break  # End of chain
            current_atom = next_atoms[0]

        # Check minimum chain length (e.g., at least 4 carbons in the chain)
        if chain_length >= 4:
            return True, "Matches 2-hydroxy fatty acid structure with sufficient aliphatic chain"

    return False, "No valid 2-hydroxy fatty acid structure found"