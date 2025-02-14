"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA is any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the CoA moiety (simplified pattern)
    coa_smarts = "NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc23)[C@H](O)[C@@H]1OP(=O)(O)O"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Find the thioester linkage (C(=O)-S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Assume the first thioester match corresponds to the fatty acyl chain
    carbonyl_c_idx = thioester_matches[0][0]
    sulfur_idx = thioester_matches[0][1]

    # Perform a BFS to extract the fatty acyl chain atoms
    acyl_chain_atoms = set()
    visited = set()
    to_visit = [carbonyl_c_idx]

    while to_visit:
        atom_idx = to_visit.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        acyl_chain_atoms.add(atom_idx)
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            nbr_idx = nbr.GetIdx()
            # Do not cross the sulfur atom into CoA moiety
            if nbr_idx == sulfur_idx:
                continue
            # Avoid entering the CoA moiety by skipping atoms beyond the sulfur
            if nbr_idx in visited:
                continue
            to_visit.append(nbr_idx)

    # Create a sub-molecule of the fatty acyl chain
    acyl_chain_mol = Chem.PathToSubmol(mol, list(acyl_chain_atoms))

    # Count the number of carbon-carbon double bonds in the acyl chain
    double_bonds = 0
    for bond in acyl_chain_mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            # Count only carbon-carbon double bonds
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                double_bonds += 1

    if double_bonds == 1:
        return True, "Contains one carbon-carbon double bond in the fatty acyl chain"
    else:
        return False, f"Contains {double_bonds} carbon-carbon double bonds in the fatty acyl chain"