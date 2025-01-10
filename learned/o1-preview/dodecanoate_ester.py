"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:<id> dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is an ester where the carboxylic acid component is lauric acid (dodecanoic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester SMARTS pattern with labels
    ester_pattern = Chem.MolFromSmarts('[#6:1](=O)[O][#6]')  # Ester group: C(=O)O-C
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check if the acyl chain is dodecanoic acid
    for match in ester_matches:
        carbonyl_c_idx = match[0]
        ester_o_idx = match[2]  # The oxygen index
        # Trace the acyl chain from the carbonyl carbon
        acyl_chain_atoms = set()
        atoms_to_visit = [carbonyl_c_idx]
        visited_atoms = set()

        while atoms_to_visit:
            atom_idx = atoms_to_visit.pop()
            if atom_idx in visited_atoms:
                continue
            visited_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            acyl_chain_atoms.add(atom_idx)

            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Skip the ester oxygen to stay on acyl chain
                if neighbor_idx == ester_o_idx:
                    continue
                bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    atoms_to_visit.append(neighbor_idx)

        # Count number of carbons in acyl chain
        carbon_count = sum(1 for idx in acyl_chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if carbon_count == 12:
            return True, "Contains lauric acid ester group"

    return False, "No lauric acid ester groups found"