"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    Any fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for fatty acyl-CoA thioester linkage
    # The fatty acyl chain is connected via a thioester linkage: C(=O)SCCN
    # We will capture the carbonyl carbon of the fatty acyl chain
    fatty_acyl_pattern = Chem.MolFromSmarts('C(=O)SCCNC(=O)')
    fatty_acyl_matches = mol.GetSubstructMatches(fatty_acyl_pattern)

    if not fatty_acyl_matches:
        return False, "No fatty acyl-CoA thioester linkage found"

    # Assume the first match is the fatty acyl chain
    fatty_acyl_match = fatty_acyl_matches[0]
    acyl_carbon_idx = fatty_acyl_match[0]  # Carbonyl carbon of the fatty acyl chain

    # Traverse the fatty acyl chain starting from the carbonyl carbon
    # We will build the acyl chain by moving from the carbonyl carbon along alkyl chain
    acyl_chain_atoms = []
    visited = set()

    def traverse_acyl_chain(atom_idx):
        """Recursively traverse the acyl chain starting from the given atom index."""
        atom = mol.GetAtomWithIdx(atom_idx)
        visited.add(atom_idx)
        if atom.GetAtomicNum() != 6:
            return  # Only consider carbon atoms

        acyl_chain_atoms.append(atom_idx)
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]

        for nbr_idx in neighbors:
            bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
            # Exclude the direction back towards the CoA moiety
            if nbr_idx in visited or nbr_idx == acyl_carbon_idx:
                continue

            # Check if the neighbor is part of the acyl chain
            # Ensure that we are not traversing into the CoA moiety
            traverse_acyl_chain(nbr_idx)

    # Start traversal from the carbon next to the carbonyl carbon
    carbonyl_carbon = mol.GetAtomWithIdx(acyl_carbon_idx)
    for neighbor in carbonyl_carbon.GetNeighbors():
        nbr_idx = neighbor.GetIdx()
        if neighbor.GetAtomicNum() == 6 and nbr_idx != fatty_acyl_match[1]:  # Exclude the sulfur atom
            # Begin traversal
            traverse_acyl_chain(nbr_idx)
            break

    # Ensure the acyl chain has at least 12 carbons
    if len(acyl_chain_atoms) < 12:
        return False, f"Fatty acyl chain has {len(acyl_chain_atoms)} carbons, less than 12"

    # Get the atom indices for C11 and C12
    # Note: acyl_chain_atoms[0] is the alpha carbon (adjacent to carbonyl carbon)
    c11_idx = acyl_chain_atoms[10]
    c12_idx = acyl_chain_atoms[11]

    # Get the bond between C11 and C12
    bond = mol.GetBondBetweenAtoms(c11_idx, c12_idx)
    if bond is None:
        return False, "No bond between C11 and C12"

    bond_type = bond.GetBondType()
    if bond_type == rdchem.BondType.SINGLE:
        return True, "C11-C12 bond is saturated (single bond)"
    elif bond_type == rdchem.BondType.DOUBLE:
        return False, "C11-C12 bond is unsaturated (double bond)"
    else:
        return False, f"C11-C12 bond is neither single nor double (bond type: {bond_type})"