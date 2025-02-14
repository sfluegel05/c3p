"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: fatty acid methyl ester
"""
from rdkit import Chem

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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the methyl ester functional group pattern
    # Pattern: Carbonyl carbon (=O) connected to an oxygen that is connected to a methyl group (ester linkage)
    methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0][CH3]")
    ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    if not ester_matches:
        return False, "No methyl ester functional group found"

    # Iterate over each methyl ester group found
    for match in ester_matches:
        carbonyl_c_idx = match[0]
        ester_o_idx = match[1]
        methyl_c_idx = match[2]

        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        ester_o = mol.GetAtomWithIdx(ester_o_idx)
        methyl_c = mol.GetAtomWithIdx(methyl_c_idx)

        # Check that the ester oxygen is connected to a methyl group (CH3)
        if methyl_c.GetDegree() != 1 or methyl_c.GetAtomicNum() != 6:
            continue  # Not a methyl group

        # Traverse the chain connected to the carbonyl carbon (excluding the ester oxygen)
        acyl_chain_atoms = set()
        atoms_to_visit = [nbr for nbr in carbonyl_c.GetNeighbors() if nbr.GetIdx() != ester_o_idx]
        while atoms_to_visit:
            atom = atoms_to_visit.pop()
            idx = atom.GetIdx()
            if idx in acyl_chain_atoms:
                continue  # Already visited

            acyl_chain_atoms.add(idx)

            # Allow common elements in fatty acid chains
            atomic_num = atom.GetAtomicNum()
            if atomic_num in [1, 6, 7, 8, 15, 16, 17, 35, 53]:
                pass  # Allowed atoms (H, C, N, O, P, S, Cl, Br, I)
            else:
                return False, f"Atom {atom.GetSymbol()} not typical in fatty acid chain"

            # Add neighbors to visit, excluding the carbonyl carbon
            neighbors = [
                nbr for nbr in atom.GetNeighbors()
                if nbr.GetIdx() != carbonyl_c_idx and nbr.GetAtomicNum() != 1
            ]
            atoms_to_visit.extend(neighbors)

        # Count the number of carbon atoms in the acyl chain
        c_count = sum(
            1 for idx in acyl_chain_atoms
            if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6
        )

        # Set a reasonable minimum chain length (e.g., 8 carbons)
        if c_count < 8:
            continue  # Chain too short to be a fatty acid

        # Passed all checks
        return True, "Contains methyl ester group with appropriate fatty acid chain"

    return False, "No suitable fatty acid methyl ester pattern found"