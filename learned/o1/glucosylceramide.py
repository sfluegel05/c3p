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
    A glucosylceramide is a ceramide (sphingoid base linked to fatty acid) with 
    a glucose moiety attached via a beta-glycosidic bond.

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

    # Check for glucose moiety (β-D-glucose)
    glucose_smarts = "[C@H]1(O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O1)CO"
    glucose_pattern = Chem.MolFromSmarts(glucose_smarts)
    if glucose_pattern is None:
        return False, "Invalid glucose SMARTS pattern"
    glucose_matches = mol.GetSubstructMatches(glucose_pattern)
    if not glucose_matches:
        return False, "No glucose moiety found"

    # Check for amide bond (ceramide backbone)
    amide_pattern = Chem.MolFromSmarts("NC(=O)C")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found"

    # Check for long-chain sphingoid base (contains chain with amino and hydroxyl groups)
    sphingoid_pattern = Chem.MolFromSmarts("[N][C@@H](CO)[C@H](O)C")
    sphingoid_matches = mol.GetSubstructMatches(sphingoid_pattern)
    if not sphingoid_matches:
        return False, "No sphingoid base found"

    # Verify that the glucose moiety is connected via a glycosidic bond
    # Find the oxygen atom connecting glucose to sphingoid base
    glycosidic_bond_found = False
    for glucose_match in glucose_matches:
        # Find the anomeric carbon (connected to two oxygens)
        anomeric_carbons = [idx for idx in glucose_match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and len([nbr for nbr in mol.GetAtomWithIdx(idx).GetNeighbors() if nbr.GetAtomicNum() == 8]) == 2]
        for anomeric_idx in anomeric_carbons:
            anomeric_atom = mol.GetAtomWithIdx(anomeric_idx)
            for neighbor in anomeric_atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if nbr_idx not in glucose_match and neighbor.GetAtomicNum() == 8:
                    # Oxygen connecting to sphingoid base
                    for nbr2 in neighbor.GetNeighbors():
                        if nbr2.GetIdx() != anomeric_idx and nbr2.GetAtomicNum() == 6:
                            # Connected carbon should be part of sphingoid base
                            sphingoid_paths = Chem.rdmolops.GetShortestPath(mol, nbr2.GetIdx(), sphingoid_matches[0][0])
                            if sphingoid_paths:
                                glycosidic_bond_found = True
                                break
                    if glycosidic_bond_found:
                        break
            if glycosidic_bond_found:
                break
        if glycosidic_bond_found:
            break
    if not glycosidic_bond_found:
        return False, "No glycosidic linkage between glucose and sphingoid base found"

    # Optional: Check that the fatty acid chain (acyl group) is long enough
    # Find the carbonyl carbon in the amide bond
    fatty_acid_length = 0
    for amide_match in amide_matches:
        carbonyl_c_idx = amide_match[2]
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        # Traverse the chain connected to the carbonyl carbon
        visited = set()
        stack = [nbr.GetIdx() for nbr in carbonyl_c.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != amide_match[1]]
        while stack:
            idx = stack.pop()
            if idx in visited:
                continue
            visited.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                fatty_acid_length += 1
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx not in visited and nbr.GetAtomicNum() == 6:
                        stack.append(nbr_idx)
    if fatty_acid_length < 12:
        return False, f"Fatty acid chain too short ({fatty_acid_length} carbons), expected at least 12 carbons"

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
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}