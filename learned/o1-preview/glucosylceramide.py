"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: CHEBI:37547 glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide is a ceramide (sphingoid base linked to a fatty acid via an amide bond) 
    with a glucose moiety attached via a beta-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glucose moiety (beta-D-glucose)
    glucose_smarts = "[C@H]1(O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O1)CO"
    glucose_pattern = Chem.MolFromSmarts(glucose_smarts)
    if not glucose_pattern:
        return False, "Invalid glucose SMARTS pattern"

    glucose_matches = mol.GetSubstructMatches(glucose_pattern)
    if not glucose_matches:
        return False, "No glucose moiety found"

    # Check for ceramide backbone (sphingoid base linked via amide bond to fatty acid)
    ceramide_smarts = "N[C@@H](CO)[C@H](O)CC=CCCCCCCCCCCCCCCCC"  # Simplified sphingosine backbone
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)
    if not ceramide_pattern:
        return False, "Invalid ceramide SMARTS pattern"
    
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if not ceramide_matches:
        return False, "No ceramide backbone found"

    # Check for amide bond connecting fatty acid to sphingoid base
    amide_pattern = Chem.MolFromSmarts("C(=O)N[C@H]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found"

    # Check for glycosidic bond between glucose and ceramide
    # Look for oxygen atom connecting glucose to sphingoid base
    glycosidic_bond_found = False
    for glucose_match in glucose_matches:
        glucose_atoms = set(glucose_match)
        # Find the anomeric carbon (connected to two oxygens)
        for idx in glucose_match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                oxy_neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
                if len(oxy_neighbors) == 2:
                    anomeric_carbon = idx
                    anomeric_oxygen = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in glucose_atoms][0]
                    # Check if this oxygen connects to ceramide backbone
                    oxygen_atom = mol.GetAtomWithIdx(anomeric_oxygen)
                    for nbr in oxygen_atom.GetNeighbors():
                        nbr_idx = nbr.GetIdx()
                        if nbr_idx not in glucose_atoms:
                            # Check if neighbor atom is part of ceramide
                            for ceramide_match in ceramide_matches:
                                if nbr_idx in ceramide_match:
                                    glycosidic_bond_found = True
                                    break
                        if glycosidic_bond_found:
                            break
                if glycosidic_bond_found:
                    break
            if glycosidic_bond_found:
                break
        if glycosidic_bond_found:
            break

    if not glycosidic_bond_found:
        return False, "No glycosidic linkage between glucose and ceramide found"

    # Optional: Check length of fatty acid chain in ceramide
    chain_lengths = []
    for amide_match in amide_matches:
        carbonyl_c_idx = amide_match[0]
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        # Find the carbon chain connected to the carbonyl carbon
        chain_length = 0
        visited = set()
        stack = []
        for nbr in carbonyl_c.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr.GetAtomicNum() == 6 and nbr_idx != amide_match[1]:
                stack.append(nbr_idx)
        while stack:
            idx = stack.pop()
            if idx in visited:
                continue
            visited.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                chain_length += 1
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx not in visited and nbr.GetAtomicNum() == 6:
                        stack.append(nbr_idx)
        chain_lengths.append(chain_length)

    if not chain_lengths or max(chain_lengths) < 12:
        return False, "Fatty acid chain too short, expected at least 12 carbons"

    return True, "Contains ceramide backbone with glucose attached via glycosidic bond"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:37547',
        'name': 'glucosylceramide',
        'definition': 'Any of the cerebrosides in which the monosaccharide head group is glucose.',
        'parents': ['CHEBI:33563']
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
    },
    'message': None,
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}