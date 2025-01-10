"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: CHEBI:XXXX fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is the carboxylic ester obtained by the formal condensation
    of a fatty acid with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the methyl ester SMARTS pattern
    methyl_ester_smarts = "[CX3](=O)[O][CH3]"  # Ester carbonyl connected to OCH3
    methyl_ester_pattern = Chem.MolFromSmarts(methyl_ester_smarts)

    # Find methyl ester groups in the molecule
    ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    if not ester_matches:
        return False, "No methyl ester group found"

    # For each methyl ester group, check if it's attached to a fatty acid chain
    for match in ester_matches:
        ester_carbon_idx = match[0]  # Carbonyl carbon atom index
        ester_oxygen_idx = match[1]  # Ester oxygen atom index (connected to methyl)
        methyl_carbon_idx = match[2]  # Methyl carbon atom index

        # Get the carbonyl carbon atom
        carbonyl_carbon = mol.GetAtomWithIdx(ester_carbon_idx)

        # Find the atom attached to the carbonyl carbon that is not the ester oxygen
        neighbors = [nbr for nbr in carbonyl_carbon.GetNeighbors() if nbr.GetIdx() != ester_oxygen_idx]
        if len(neighbors) != 1:
            continue  # Skip if not exactly one neighbor (should be the fatty acid chain)
        chain_start_atom = neighbors[0]

        # Traverse the chain to check if it's a fatty acid chain
        chain_atom_indices = set()
        atoms_to_visit = [chain_start_atom.GetIdx()]
        while atoms_to_visit:
            current_idx = atoms_to_visit.pop()
            if current_idx in chain_atom_indices:
                continue
            chain_atom_indices.add(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)
            atomic_num = current_atom.GetAtomicNum()

            # Allow carbon, hydrogen, and oxygen atoms (to include common fatty acid functional groups)
            if atomic_num not in [1, 6, 8]:  # H, C, O
                # Disallow other heteroatoms
                return False, f"Chain contains disallowed atom with atomic number {atomic_num}"

            # Allow small rings (e.g., epoxides) by not excluding ring atoms
            # Add neighbor atoms to visit
            for nbr in current_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx != ester_carbon.GetIdx() and nbr_idx not in chain_atom_indices:
                    atoms_to_visit.append(nbr_idx)

        # Estimate chain length by counting carbon atoms in the chain
        chain_carbons = [idx for idx in chain_atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        if len(chain_carbons) < 4:
            # Too short to be a fatty acid
            continue

        # If we reach here, we have found a fatty acid methyl ester
        return True, "Molecule is a fatty acid methyl ester"

    # No suitable fatty acid methyl ester found
    return False, "No suitable fatty acid methyl ester group found"