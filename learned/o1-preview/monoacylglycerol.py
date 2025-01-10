"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: CHEBI:35700 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol is a glyceride in which any one of the hydroxyl groups is esterified with an acyl group,
    and the remaining two hydroxyl groups can be either unmodified (OH) or alkyl ethers (OR).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphorus atoms (exclude phospholipids)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Contains phosphorus atom, not a monoacylglycerol"

    # Check for sulfur atoms (exclude sulfolipids)
    if any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms()):
        return False, "Contains sulfur atom, not a monoacylglycerol"

    # Identify ester bonds: [#6](=O)[O][#6]
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Contains {len(ester_matches)} ester bonds, expected exactly 1"

    ester_bond_atoms = ester_matches[0]
    ester_oxygen_idx = ester_bond_atoms[2]

    # Identify glycerol backbone: [CH2][CH][CH2]
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # Check that the ester bond is connected to the glycerol backbone
    ester_oxygen = mol.GetAtomWithIdx(ester_oxygen_idx)

    connected_to_glycerol = False
    for glycerol_match in glycerol_matches:
        glycerol_atoms = [mol.GetAtomWithIdx(idx) for idx in glycerol_match]
        glycerol_atom_indices = glycerol_match

        # Check if ester oxygen is connected to one of the glycerol carbons
        for atom in glycerol_atoms:
            if ester_oxygen.IsInNeighbor(atom):
                connected_to_glycerol = True
                glycerol_carbon_with_ester = atom.GetIdx()

                # Check other two carbons have hydroxyl or ether groups
                other_carbons = [idx for idx in glycerol_atom_indices if idx != glycerol_carbon_with_ester]
                hydroxyl_or_ether_count = 0
                for idx in other_carbons:
                    carbon = mol.GetAtomWithIdx(idx)
                    has_oxygen = False
                    for nbr in carbon.GetNeighbors():
                        if nbr.GetAtomicNum() == 8:
                            oxygen = nbr
                            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), oxygen.GetIdx())
                            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                has_oxygen = True
                    if not has_oxygen:
                        break  # Missing hydroxyl or ether group
                    else:
                        hydroxyl_or_ether_count +=1
                if hydroxyl_or_ether_count != 2:
                    continue  # Try next glycerol match
                else:
                    return True, "Contains glycerol backbone with one esterified hydroxyl group and two hydroxyl/ether groups"

    if not connected_to_glycerol:
        return False, "Ester bond not connected to glycerol backbone or glycerol backbone missing hydroxyl groups"

    return False, "Does not meet criteria for monoacylglycerol"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35700',
        'name': 'monoacylglycerol',
        'definition': 'A glyceride in which any one of the R groups (position not specified) is an acyl group while the remaining two R groups can be either H or alkyl groups.',
        'parents': ['CHEBI:25212', 'CHEBI:17855']
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