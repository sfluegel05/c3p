"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    A sphingomyelin is defined as a phospholipid where the amino group of a sphingoid base is in amide linkage with one of several fatty acids,
    while the terminal hydroxy group of the sphingoid base is esterified to phosphorylcholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphorylcholine group attached via ester linkage
    phosphocholine_pattern = Chem.MolFromSmarts("O[P](=O)(OCC[N+](C)(C)C)O")
    phospho_matches = mol.GetSubstructMatches(phosphocholine_pattern)
    if not phospho_matches:
        return False, "No phosphorylcholine group found"

    # Check for amide linkage with fatty acid
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide linkage found between amino group and fatty acid"

    # For each amide linkage, verify sphingoid base and fatty acid
    for amide_match in amide_matches:
        carbonyl_c_idx = amide_match[0]
        amide_n_idx = amide_match[2]
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        amide_n = mol.GetAtomWithIdx(amide_n_idx)

        # Check fatty acid chain length (carbons connected to carbonyl carbon)
        fatty_acid_length = get_chain_length(carbonyl_c, exclude_idxs={amide_n_idx})
        if fatty_acid_length < 12:
            continue  # Fatty acid chain too short

        # Check sphingoid base chain length (carbons connected to amide nitrogen)
        sphingoid_chain_length = get_chain_length(amide_n, exclude_idxs={carbonyl_c_idx})
        if sphingoid_chain_length < 12:
            continue  # Sphingoid base chain too short

        # Check for hydroxyl group(s) on sphingoid base
        hydroxyls_on_sphingoid = count_hydroxyls(amide_n, exclude_idxs={carbonyl_c_idx})
        if hydroxyls_on_sphingoid < 1:
            continue  # No hydroxyl groups on sphingoid base

        # Check for phosphorylcholine attached to sphingoid base hydroxyl
        phospho_attached = False
        for match in phospho_matches:
            phospho_o_idx = match[0]  # Oxygen connected to phosphate
            phospho_o = mol.GetAtomWithIdx(phospho_o_idx)
            # Check if this oxygen is connected to sphingoid base
            if is_connected(phospho_o, amide_n, exclude_idxs={carbonyl_c_idx}):
                phospho_attached = True
                break
        if not phospho_attached:
            continue  # Phosphorylcholine not attached to sphingoid base

        # All checks passed
        return True, "Molecule is a sphingomyelin with sphingoid base, amide-linked fatty acid, and phosphorylcholine group"

    # If none of the amide linkages satisfy the conditions
    return False, "No sphingomyelin structure found"

def get_chain_length(start_atom, exclude_idxs):
    """
    Counts the number of carbon atoms in the chain starting from the given atom.

    Args:
        start_atom (Atom): Starting RDKit atom.
        exclude_idxs (set): Atom indices to exclude from traversal.

    Returns:
        int: Number of carbon atoms in the chain.
    """
    visited = set()
    to_visit = [start_atom]
    chain_length = 0
    while to_visit:
        atom = to_visit.pop()
        idx = atom.GetIdx()
        if idx in visited or idx in exclude_idxs:
            continue
        visited.add(idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_length += 1
            neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in visited]
            to_visit.extend(neighbors)
    return chain_length

def count_hydroxyls(start_atom, exclude_idxs):
    """
    Counts the number of hydroxyl groups connected to the chain starting from the given atom.

    Args:
        start_atom (Atom): Starting RDKit atom.
        exclude_idxs (set): Atom indices to exclude from traversal.

    Returns:
        int: Number of hydroxyl groups.
    """
    visited = set()
    to_visit = [start_atom]
    hydroxyl_count = 0
    while to_visit:
        atom = to_visit.pop()
        idx = atom.GetIdx()
        if idx in visited or idx in exclude_idxs:
            continue
        visited.add(idx)
        # Check for hydroxyl group (oxygen with single bond to hydrogen)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # Oxygen
                for onbr in nbr.GetNeighbors():
                    if onbr.GetAtomicNum() == 1:  # Hydrogen
                        hydroxyl_count += 1
        neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in visited]
        to_visit.extend(neighbors)
    return hydroxyl_count

def is_connected(atom1, atom2, exclude_idxs):
    """
    Determines if two atoms are connected within the molecule.

    Args:
        atom1 (Atom): First RDKit atom.
        atom2 (Atom): Second RDKit atom.
        exclude_idxs (set): Atom indices to exclude from traversal.

    Returns:
        bool: True if atoms are connected, False otherwise.
    """
    visited = set()
    to_visit = [atom1]
    while to_visit:
        atom = to_visit.pop()
        idx = atom.GetIdx()
        if idx in visited or idx in exclude_idxs:
            continue
        if idx == atom2.GetIdx():
            return True
        visited.add(idx)
        neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in visited]
        to_visit.extend(neighbors)
    return False