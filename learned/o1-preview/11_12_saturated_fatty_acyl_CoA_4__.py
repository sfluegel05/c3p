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
    This includes all fatty acyl-CoA(4-) molecules where the bond between the 11th and 12th carbons
    in the fatty acyl chain is either saturated (single bond) or absent (for chains shorter than 12 carbons).

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
    # The fatty acyl chain is connected via a thioester linkage: C(=O)SCCNC(=O)
    # We will capture the carbonyl carbon of the fatty acyl chain
    fatty_acyl_pattern = Chem.MolFromSmarts('C(=O)SCCNC(=O)')
    fatty_acyl_matches = mol.GetSubstructMatches(fatty_acyl_pattern)

    if not fatty_acyl_matches:
        return False, "No fatty acyl-CoA thioester linkage found"

    # Assume the first match is the fatty acyl chain
    fatty_acyl_match = fatty_acyl_matches[0]
    acyl_carbon_idx = fatty_acyl_match[0]  # Carbonyl carbon of the fatty acyl chain

    # Get the acyl chain starting from the alpha carbon
    def get_acyl_chain(mol, acyl_carbon_idx):
        """Get the acyl chain starting from the alpha carbon."""
        # Get the alpha carbon(s)
        carbonyl_carbon = mol.GetAtomWithIdx(acyl_carbon_idx)
        alpha_carbons = [nbr for nbr in carbonyl_carbon.GetNeighbors() 
                         if nbr.GetAtomicNum() == 6 and
                         mol.GetBondBetweenAtoms(acyl_carbon_idx, nbr.GetIdx()).GetBondType() == rdchem.BondType.SINGLE]
        if not alpha_carbons:
            return [], "No alpha carbon found"
        if len(alpha_carbons) > 1:
            return [], "Multiple alpha carbons found, branching detected"

        alpha_carbon_idx = alpha_carbons[0].GetIdx()
        chain = [acyl_carbon_idx, alpha_carbon_idx]

        previous_atom_idx = acyl_carbon_idx
        current_atom_idx = alpha_carbon_idx

        while True:
            current_atom = mol.GetAtomWithIdx(current_atom_idx)
            # Get neighboring carbons excluding the previous atom
            neighbors = [nbr for nbr in current_atom.GetNeighbors() 
                         if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != previous_atom_idx]
            if not neighbors:
                break  # Reached end of chain
            if len(neighbors) > 1:
                return [], f"Branching detected at atom {current_atom_idx}"
            next_atom_idx = neighbors[0].GetIdx()
            chain.append(next_atom_idx)
            previous_atom_idx = current_atom_idx
            current_atom_idx = next_atom_idx

        return chain, None

    acyl_chain, error = get_acyl_chain(mol, acyl_carbon_idx)
    if error:
        return False, error

    # Count the number of carbons in the acyl chain (including carbonyl carbon)
    num_carbons = len(acyl_chain)
    # For chains shorter than 12 carbons, we consider the C11-C12 bond as saturated (since it does not exist)
    if num_carbons < 13:
        return True, f"Fatty acyl chain has {num_carbons - 1} carbons, C11-C12 bond is absent and considered saturated"

    # Get the atom indices for C11 and C12 (0-based indexing, acyl_carbon_idx is C1)
    c11_idx = acyl_chain[11]
    c12_idx = acyl_chain[12]

    # Get the bond between C11 and C12
    bond = mol.GetBondBetweenAtoms(c11_idx, c12_idx)
    if bond is None:
        return False, "No bond between C11 and C12"

    bond_type = bond.GetBondType()
    if bond_type == rdchem.BondType.SINGLE:
        return True, "C11-C12 bond is saturated (single bond)"
    elif bond_type == rdchem.BondType.DOUBLE or bond_type == rdchem.BondType.AROMATIC:
        return False, "C11-C12 bond is unsaturated (double bond)"
    else:
        return False, f"C11-C12 bond is neither single nor double bond (bond type: {bond_type})"