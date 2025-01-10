"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: fatty amide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define amide functional group pattern: [C](=O)N
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide functional group found"

    # For each amide group, check the acyl chain attached to carbonyl carbon
    for match in amide_matches:
        carbonyl_c_index = match[0]
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_index)

        # Get the chain attached to the carbonyl carbon (excluding the carbonyl oxygen)
        acyl_chain_atoms = set()
        atoms_to_visit = [neighbor for neighbor in carbonyl_c.GetNeighbors() if neighbor.GetAtomicNum() == 6]
        visited_atoms = set()
        while atoms_to_visit:
            atom = atoms_to_visit.pop()
            atom_idx = atom.GetIdx()
            if atom_idx in visited_atoms:
                continue
            visited_atoms.add(atom_idx)
            acyl_chain_atoms.add(atom_idx)

            # Exclude atoms in rings
            if atom.IsInRing():
                continue

            # Exclude heteroatoms
            if atom.GetAtomicNum() != 6:
                continue

            # Add neighboring carbons to visit
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited_atoms and neighbor.GetAtomicNum() == 6:
                    atoms_to_visit.append(neighbor)

        # Count the number of carbons in the acyl chain
        acyl_chain_length = len(acyl_chain_atoms)
        if acyl_chain_length >= 4:
            return True, "Molecule is a fatty amide with acyl chain derived from a fatty acid"

    return False, "No suitable acyl chain found for fatty amide classification"


__metadata__ = {   'chemical_class': {   'name': 'fatty amide',
                                  'definition': 'A monocarboxylic acid amide derived from a fatty acid.'},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 1,
    'success': True}