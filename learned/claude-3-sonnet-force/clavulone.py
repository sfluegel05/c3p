"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: CHEBI:137276 clavulone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondStereo

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    Clavulones are esterified prostanoids obtained from marine corals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclopentenone core
    cyclopentenone_pattern = Chem.MolFromSmarts("[C@@]1(=O)CC=CC1")
    if not mol.HasSubstructMatch(cyclopentenone_pattern):
        return False, "No cyclopentenone core found"

    # Look for trans double bonds in side chains
    trans_double_bond_pattern = Chem.MolFromSmarts("/C=C/CCCC")
    trans_double_bond_matches = mol.GetSubstructMatches(trans_double_bond_pattern)
    if not trans_double_bond_matches:
        return False, "No trans double bonds in side chains"

    # Look for halogen (Cl, Br, I) attached to the cyclopentenone ring
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I]")
    halogen_matches = mol.GetSubstructMatches(halogen_pattern, useQueryQueryMatches=True)
    halogen_on_ring = False
    for match in halogen_matches:
        atom_idx = match[0]
        if mol.GetAtomWithIdx(atom_idx).IsInRing():
            halogen_on_ring = True
            break
    if not halogen_on_ring:
        return False, "No halogen substituent on the cyclopentenone ring"

    # Look for epoxide ring fused to the cyclopentenone ring
    epoxide_pattern = Chem.MolFromSmarts("[C@]12OC[C@@]1(O)CC2=O")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Check for specific positions and stereochemistry of ester groups
    # (Attached to C10 and C12 of prostanoid backbone)
    ester_positions = []
    for match in ester_matches:
        atom_idx = match[1]  # Ester carbon atom index
        atom = mol.GetAtomWithIdx(atom_idx)
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                neighbor_atom = bond.GetOtherAtom(atom)
                if neighbor_atom.GetDegree() == 4 and neighbor_atom.GetTotalNumHs() == 1:
                    ester_positions.append(neighbor_atom.GetIdx())
    expected_ester_positions = [10, 12]  # Assuming standard prostanoid numbering
    if set(ester_positions) != set(expected_ester_positions):
        return False, "Incorrect positions or stereochemistry of ester groups"

    # Check for long carbon chains (>5 carbons)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if not long_chain_matches:
        return False, "No long carbon chains found"

    # Additional checks based on molecular properties
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 800:
        return False, "Molecular weight outside typical range for clavulones"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 20 or o_count < 5:
        return False, "Insufficient carbon or oxygen atoms for clavulones"

    # If all checks pass, classify as clavulone
    if epoxide_matches:
        return True, "Contains cyclopentenone core, trans double bonds, halogen on ring, epoxide ring, and ester groups in correct positions"
    else:
        return True, "Contains cyclopentenone core, trans double bonds, halogen on ring, and ester groups in correct positions"

__metadata__ = {'chemical_class': {'id': 'CHEBI:137276',
                                   'name': 'clavulone',
                                   'definition': 'A class of esterified prostanoids obtained from marine corals.',
                                   'parents': ['CHEBI:23058', 'CHEBI:38513', 'CHEBI:50114']},
                 'config': {'llm_model_name': 'lbl/claude-sonnet',
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
                 'success': True,
                 'best': True,
                 'error': '',
                 'stdout': None,
                 'num_true_positives': 184,
                 'num_false_positives': 0,
                 'num_true_negatives': 182400,
                 'num_false_negatives': 0,
                 'num_negatives': None,
                 'precision': 1.0,
                 'recall': 1.0,
                 'f1': 1.0,
                 'accuracy': 1.0}