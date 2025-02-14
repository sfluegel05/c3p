"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester is a fatty acid ester resulting from the formal condensation
    of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester functional group SMARTS pattern
    ester_smarts = '[C;!R](=O)O[C;!R]'  # Non-ring ester functional group
    ester_mol = Chem.MolFromSmarts(ester_smarts)

    ester_matches = mol.GetSubstructMatches(ester_mol)
    if not ester_matches:
        return False, "No ester groups found"

    # Iterate over each ester group found
    for match in ester_matches:
        carbonyl_carbon_idx = match[0]
        ester_oxygen_idx = match[2]

        # Traverse the acyl side (decanoic acid side)
        acyl_atoms = set()
        atoms_to_visit = [carbonyl_carbon_idx]
        while atoms_to_visit:
            current_atom_idx = atoms_to_visit.pop()
            if current_atom_idx in acyl_atoms:
                continue
            acyl_atoms.add(current_atom_idx)
            current_atom = mol.GetAtomWithIdx(current_atom_idx)
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                bond = mol.GetBondBetweenAtoms(current_atom_idx, neighbor_idx)
                if bond.GetBondType() == rdchem.BondType.SINGLE and neighbor_idx != ester_oxygen_idx:
                    atoms_to_visit.append(neighbor_idx)

        # Count the number of carbon atoms in the acyl chain
        carbon_count = sum(1 for idx in acyl_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)

        # Check if the acyl chain corresponds to decanoic acid (10 carbons)
        if carbon_count == 10:
            return True, "Contains decanoate ester group"

    return False, "No decanoate ester groups with 10 carbons found"