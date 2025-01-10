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

    # Exclude molecules with multiple amide groups (e.g., peptides)
    if len(amide_matches) > 1:
        return False, "Multiple amide groups found, may be peptide"

    # For each amide group, check the acyl chain length
    for match in amide_matches:
        carbonyl_c_index = match[0]
        nitrogen_index = match[2]

        # Get the carbon chain attached to the carbonyl carbon (acyl chain)
        acyl_chain_atoms = []
        visited_atoms = set()
        atoms_to_visit = [mol.GetAtomWithIdx(carbonyl_c_index)]

        while atoms_to_visit:
            atom = atoms_to_visit.pop()
            if atom.GetIdx() in visited_atoms:
                continue
            visited_atoms.add(atom.GetIdx())
            atomic_num = atom.GetAtomicNum()
            # Skip the carbonyl carbon and heteroatoms
            if atom.GetIdx() != carbonyl_c_index and atomic_num != 6:
                continue
            acyl_chain_atoms.append(atom)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != carbonyl_c_index:
                    atoms_to_visit.append(neighbor)

        # Count the number of carbons in the acyl chain
        acyl_chain_length = len(acyl_chain_atoms)
        if acyl_chain_length < 8:
            return False, f"Acyl chain too short ({acyl_chain_length} carbons), not a fatty amide"

        # Check for branching in the acyl chain
        for atom in acyl_chain_atoms:
            if atom.GetDegree() > 2 and atom.GetIdx() != carbonyl_c_index:
                return False, "Branched acyl chain, not a fatty amide"

    return True, "Molecule is a fatty amide with a long acyl chain derived from a fatty acid"


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
    'attempt': 0,
    'success': True}