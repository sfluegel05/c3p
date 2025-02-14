"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: CHEBI:76955 N-acylphytosphingosine
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide that is phytosphingosine
    having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for phytosphingosine backbone
    # Phytosphingosine backbone: C-C-N-C-C with hydroxyl groups at positions 1, 3, and 4
    phytosphingosine_smarts = "[C@@H](O)[C@H](O)[C@H](CO)N"
    phytosphingosine_pattern = Chem.MolFromSmarts(phytosphingosine_smarts)
    if not phytosphingosine_pattern:
        return False, "Failed to create phytosphingosine SMARTS pattern"

    # Check for phytosphingosine backbone
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"

    # Define SMARTS pattern for N-acyl group (amide bond)
    n_acyl_smarts = "N[C](=O)[C;H2,C;H3]"
    n_acyl_pattern = Chem.MolFromSmarts(n_acyl_smarts)
    if not n_acyl_pattern:
        return False, "Failed to create N-acyl SMARTS pattern"

    # Check for N-acylation
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl group attached to nitrogen"

    # Verify long aliphatic chains
    # Count the number of carbons in the fatty acyl chain and sphingoid base chain

    # Find the amide bond
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found"

    # Assume first amide bond is the N-acylation site
    carbonyl_c_idx, nitrogen_idx = amide_matches[0]

    # Trace the fatty acyl chain from the carbonyl carbon
    fatty_acyl_length = 0
    visited_atoms = set()
    atoms_to_visit = [carbonyl_c_idx]
    while atoms_to_visit:
        atom_idx = atoms_to_visit.pop()
        if atom_idx in visited_atoms:
            continue
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon
            fatty_acyl_length += 1
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() not in visited_atoms and nbr.GetIdx() != nitrogen_idx]
            atoms_to_visit.extend(neighbors)

    # Check if fatty acyl chain is sufficiently long (e.g., at least 10 carbons)
    if fatty_acyl_length < 10:
        return False, f"Fatty acyl chain is too short ({fatty_acyl_length} carbons)"

    # Trace the sphingoid base chain from the alpha carbon (next to nitrogen)
    sphingoid_chain_length = 0
    visited_atoms = set()
    alpha_carbon_idx = None
    nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
    for neighbor in nitrogen_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != carbonyl_c_idx:
            alpha_carbon_idx = neighbor.GetIdx()
            break
    if alpha_carbon_idx is None:
        return False, "Cannot find alpha carbon next to nitrogen"

    atoms_to_visit = [alpha_carbon_idx]
    while atoms_to_visit:
        atom_idx = atoms_to_visit.pop()
        if atom_idx in visited_atoms:
            continue
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon
            sphingoid_chain_length += 1
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() not in visited_atoms and nbr.GetAtomicNum() == 6]
            atoms_to_visit.extend(neighbors)

    # Check if sphingoid base chain is sufficiently long (e.g., at least 15 carbons)
    if sphingoid_chain_length < 15:
        return False, f"Sphingoid base chain is too short ({sphingoid_chain_length} carbons)"

    # Check for hydroxyl groups at correct positions
    # Positions 1, 3, and 4 relative to sphingoid base
    hydroxyl_counts = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    if mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                        hydroxyl_counts += 1
    if hydroxyl_counts < 3:
        return False, f"Not enough hydroxyl groups ({hydroxyl_counts} found)"

    return True, "Molecule matches N-acylphytosphingosine pattern"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:76955',
        'name': 'N-acylphytosphingosine',
        'definition': 'A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.',
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    }
}