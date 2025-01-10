"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: CHEBI:12348 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is a molecule consisting of a cyclohexane ring in the myo-configuration,
    with hydroxyl and phosphate groups attached to the ring carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find cyclohexane rings (6-membered carbon rings)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    cyclohexane_ring = None
    for ring in atom_rings:
        if len(ring) == 6:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            if all(atom.GetAtomicNum() == 6 for atom in atoms_in_ring):
                cyclohexane_ring = ring
                break
    if cyclohexane_ring is None:
        return False, "No cyclohexane ring of size 6 found"

    ring_atoms = cyclohexane_ring

    # Prepare SMARTS patterns
    hydroxyl_smarts = Chem.MolFromSmarts('[CX4;H1][OX2H]')  # Carbon with single bond to hydroxyl
    phosphate_smarts = Chem.MolFromSmarts('[CX4;H1][O][P](=O)([O])[O]')  # Carbon with single bond to phosphate group

    # Initialize flags
    phosphate_found = False

    # Check substituents on ring carbons
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        non_ring_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring_atoms]
        if not non_ring_neighbors:
            return False, "Ring carbon has no substituents"

        substituent_types = []
        has_valid_substituent = False
        for nbr in non_ring_neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # Create a fragment for matching
            atom_idxs = [atom.GetIdx(), nbr.GetIdx()]
            # Include atoms up to two bonds away for matching phosphate groups
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetIdx() != atom.GetIdx():
                    atom_idxs.append(nbr2.GetIdx())
                    for nbr3 in nbr2.GetNeighbors():
                        if nbr3.GetIdx() not in atom_idxs:
                            atom_idxs.append(nbr3.GetIdx())
            frag = Chem.PathToSubmol(mol, atom_idxs)
            if frag.HasSubstructMatch(hydroxyl_smarts):
                substituent_types.append('hydroxyl')
                has_valid_substituent = True
            elif frag.HasSubstructMatch(phosphate_smarts):
                substituent_types.append('phosphate')
                phosphate_found = True
                has_valid_substituent = True
            else:
                continue  # Could be hydrogen or other small atom

        if not has_valid_substituent:
            return False, "Ring carbon has substituents other than hydroxyl or phosphate groups"

    if not phosphate_found:
        return False, "No phosphate group found in molecule"

    # Check for extra atoms (allow for counterions like sodium)
    mol_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    allowed_elements = {'C', 'H', 'O', 'P'}
    if not all(elem in allowed_elements for elem in mol_formula if elem.isalpha()):
        return False, "Molecule contains elements other than C, H, O, P"

    return True, "Molecule matches myo-inositol phosphate structure"

# Metadata
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:12348',
        'name': 'myo-inositol phosphate',
        'definition': 'An inositol phosphate in which the inositol component has myo-configuration.',
        'parents': ['CHEBI:24848', 'CHEBI:24845']
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