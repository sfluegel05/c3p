"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: medium-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any medium-chain fatty acid.
"""

from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the thioester linkage to CoA: C(=O)SCCNC(=O)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCN")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) < 1:
        return False, "Thioester linkage to CoA moiety not found"

    # Get the carbonyl carbon atom index in the thioester linkage
    carbonyl_c_idx = thioester_matches[0][0]
    sulfur_idx = thioester_matches[0][2]

    # Traverse the acyl chain starting from the carbonyl carbon, excluding the path towards CoA
    acyl_chain_atoms = set()
    atoms_to_visit = [(carbonyl_c_idx, sulfur_idx)]
    visited_atoms = set()

    while atoms_to_visit:
        current_atom_idx, parent_atom_idx = atoms_to_visit.pop()
        if current_atom_idx in visited_atoms:
            continue
        visited_atoms.add(current_atom_idx)
        atom = mol.GetAtomWithIdx(current_atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            # Check if atom is aromatic; skip if it is
            if atom.GetIsAromatic():
                return False, "Aromatic atoms found in acyl chain, not a fatty acid"
            acyl_chain_atoms.add(current_atom_idx)
        # Traverse neighbors
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx != parent_atom_idx:
                neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
                # Exclude heteroatoms (non-carbon atoms) and atoms towards CoA (sulfur atom)
                if neighbor_atom.GetAtomicNum() == 6 and not neighbor_atom.GetIsAromatic():
                    bond = mol.GetBondBetweenAtoms(current_atom_idx, neighbor_idx)
                    # Only traverse single or double bonds (exclude triple bonds)
                    if bond.GetBondType() in [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE]:
                        atoms_to_visit.append((neighbor_idx, current_atom_idx))
                else:
                    continue  # Stop traversal if a non-carbon atom is encountered

    # Count the number of carbons in the acyl chain
    num_carbons = len(acyl_chain_atoms)

    # Check if the number of carbons is within medium-chain length (6-14)
    if num_carbons < 6:
        return False, f"Acyl chain too short ({num_carbons} carbons), not medium-chain"
    if num_carbons > 14:
        return False, f"Acyl chain too long ({num_carbons} carbons), not medium-chain"

    # Verify the presence of coenzyme A moiety
    # Define a SMARTS pattern for key features of CoA
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](N2C=NC3=C(N)N=CN=C23)[C@H](OP(=O)(O))[C@@H]1O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"

    return True, f"Contains medium-chain acyl group ({num_carbons} carbons) attached to coenzyme A via thioester linkage"