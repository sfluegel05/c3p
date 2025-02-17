"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: Fatty acid methyl ester
A fatty acid methyl ester is defined as a fatty acid ester that is the carboxylic ester obtained
by the formal condensation of a fatty acid with methanol.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    The approach is to search for a methyl ester group [CX3](=O)O[CH3] and then verify that
    the acyl portion (the fatty acid part) is reasonably long (by counting at least 3 connected carbon atoms).

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

    # Define SMARTS pattern for a methyl ester group: carbonyl carbon bonded to an O bound to CH3.
    # This represents the substructure C(=O)OCH3
    ester_smarts = "[CX3](=O)O[CH3]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return False, "Error defining ester SMARTS"

    matches = mol.GetSubstructMatches(ester_pattern)
    if not matches:
        return False, "No methyl ester group found in the molecule"

    # Helper function: count number of carbon atoms in a connected substructure.
    def count_connected_carbons(start_atom, exclude_idx):
        """
        Counts carbon atoms (atomic number 6) in the connected subgraph starting from start_atom.
        We ignore the atom with index exclude_idx to avoid going back into the ester group.
        """
        count = 0
        visited = set()
        queue = [start_atom]
        while queue:
            atom = queue.pop(0)
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6:
                count += 1
                # Enqueue neighbors that are carbons.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() == exclude_idx:
                        continue
                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                        queue.append(nbr)
        return count

    # For each methyl ester match, validate the structure further.
    # Each match is a tuple of atom indices corresponding to the SMARTS pattern:
    # (carbonyl carbon, carbonyl oxygen, ester oxygen, methyl carbon)
    for match in matches:
        carbonyl_idx, carbonylO_idx, esterO_idx, methyl_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Identify the acyl side: the neighbor of carbonyl_atom that is not the carbonyl oxygen.
        acyl_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == carbonylO_idx:
                continue
            # We expect the fatty acid chain to be attached via a carbon atom.
            if nbr.GetAtomicNum() == 6:
                acyl_neighbor = nbr
                break
        if acyl_neighbor is None:
            continue  # This match does not have an acyl chain; try next.

        # Count the number of carbon atoms in the connected acyl fragment.
        # We exclude the carbonyl carbon (to avoid counting it twice).
        chain_length = count_connected_carbons(acyl_neighbor, carbonyl_idx)

        # We use a threshold of at least 3 carbon atoms in the acyl portion.
        if chain_length < 3:
            # The acyl chain appears too short to be considered a fatty acid.
            continue

        # Additional verification: we can inspect a couple of additional molecular descriptors.
        # For example, many fatty acid methyl esters show a moderate number of rotatable bonds.
        # Also, we expect the methyl ester group to be terminal.
        n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        
        reason = (f"Found a methyl ester group with an acyl chain of {chain_length} carbon(s); "
                  f"rotatable bonds: {n_rotatable}, molecular weight: {mol_wt:.1f} Da")
        return True, reason

    return False, "Molecule contains a methyl ester group but lacks a sufficiently long acyl (fatty acid) chain"