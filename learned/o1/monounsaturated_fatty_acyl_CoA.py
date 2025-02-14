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

    # Define a SMARTS pattern for the CoA moiety (simplified but more accurate)
    coa_smarts = """
    N([H])C(=O)C(C)N([H])C(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O
    [C@H](N2C=NC3=C(N)N=CN=C32)[C@H](O)[C@@H]1OP(=O)(O)O
    """
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Invalid CoA SMARTS pattern"

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Find the thioester linkage (C(=O)-S-C)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S[C!#16]")  # Exclude disulfide bonds
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Assume the thioester linkage corresponds to the fatty acyl chain
    for match in thioester_matches:
        carbonyl_c_idx = match[0]
        sulfur_idx = match[1]
        acyl_c_idx = match[2]

        # Verify that the sulfur atom is connected to the CoA moiety
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        coa_matches = mol.GetSubstructMatches(coa_pattern)
        coa_atoms = set()
        for coa_match in coa_matches:
            coa_atoms.update(coa_match)
        if sulfur_idx not in coa_atoms:
            continue  # This sulfur is not part of CoA, skip

        # Perform a traversal from the carbonyl carbon to extract the fatty acyl chain
        acyl_chain_atoms = set()
        visited = set()
        to_visit = [carbonyl_c_idx]

        while to_visit:
            atom_idx = to_visit.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)

            # Stop traversal if we reach the sulfur atom or CoA moiety
            if atom_idx == sulfur_idx or atom_idx in coa_atoms:
                continue

            acyl_chain_atoms.add(atom_idx)

            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if nbr_idx not in visited:
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

        # Check if the acyl chain is predominantly composed of carbons (fatty acids are long hydrocarbon chains)
        c_count = sum(1 for atom in acyl_chain_mol.GetAtoms() if atom.GetAtomicNum() == 6)
        h_count = sum(1 for atom in acyl_chain_mol.GetAtoms() if atom.GetAtomicNum() == 1)
        heteroatom_count = len(acyl_chain_mol.GetAtoms()) - c_count - h_count

        if heteroatom_count > 0:
            return False, "Acyl chain contains heteroatoms, not a fatty acyl chain"

        # Decide based on the number of double bonds
        if double_bonds == 1:
            return True, "Contains one carbon-carbon double bond in the fatty acyl chain"
        elif double_bonds == 0:
            return False, "No carbon-carbon double bonds in the fatty acyl chain"
        else:
            return False, f"Contains {double_bonds} carbon-carbon double bonds in the fatty acyl chain"

    return False, "No suitable thioester linkage associated with CoA found"