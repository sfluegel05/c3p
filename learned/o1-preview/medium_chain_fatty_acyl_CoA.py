"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: medium-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any medium-chain fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

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

    # Identify the thioester linkage to CoA: C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) == 0:
        return False, "Thioester linkage to CoA not found"

    # Assume the first thioester linkage is the fatty acyl attachment point
    carbonyl_c_idx = thioester_matches[0][0]
    sulfur_idx = thioester_matches[0][2]

    # Traverse the acyl chain starting from the carbonyl carbon
    acyl_chain_atoms = set()
    atoms_to_visit = [carbonyl_c_idx]
    visited_atoms = set()

    while atoms_to_visit:
        current_atom_idx = atoms_to_visit.pop()
        if current_atom_idx in visited_atoms:
            continue
        visited_atoms.add(current_atom_idx)
        atom = mol.GetAtomWithIdx(current_atom_idx)
        atom_num = atom.GetAtomicNum()
        if atom_num == 6:  # Carbon atom
            acyl_chain_atoms.add(current_atom_idx)
        elif atom_num != 1:  # Stop if not carbon or hydrogen
            continue
        # Traverse neighbors excluding the path towards the sulfur atom
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx != sulfur_idx and neighbor_idx not in visited_atoms:
                bond = mol.GetBondBetweenAtoms(current_atom_idx, neighbor_idx)
                # Only traverse single or double bonds
                if bond.GetBondType() in [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE]:
                    atoms_to_visit.append(neighbor_idx)

    # Count the number of carbons in the acyl chain
    num_carbons = len(acyl_chain_atoms)
    # Adjust for medium-chain length (6-12 carbons)
    if num_carbons < 6:
        return False, f"Acyl chain too short ({num_carbons} carbons), not medium-chain"
    if num_carbons > 12:
        return False, f"Acyl chain too long ({num_carbons} carbons), not medium-chain"

    # Verify the presence of Coenzyme A moiety by checking for key substructures

    # Adenine ring pattern
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)N")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Adenine moiety not found, not CoA"

    # Ribose sugar connected to adenine
    ribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](CO)O1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "Ribose sugar moiety not found, not CoA"

    # Diphosphate linkage (O-P(=O)(O)-O-P(=O)(O)-O)
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "Diphosphate linkage not found, not CoA"

    # Pantetheine moiety pattern (NC(=O)CCNC(=O)C)
    pantetheine_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)C[C@H](O)C(C)(C)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Pantetheine moiety not found, not CoA"

    # Confirm that the sulfur atom is part of the pantetheine moiety
    sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
    sulfur_neighbors = [n.GetAtomicNum() for n in sulfur_atom.GetNeighbors()]
    if 7 not in sulfur_neighbors:
        return False, "Sulfur atom not connected to pantetheine moiety, not CoA"

    return True, f"Contains medium-chain acyl group ({num_carbons} carbons) attached to Coenzyme A via thioester linkage"