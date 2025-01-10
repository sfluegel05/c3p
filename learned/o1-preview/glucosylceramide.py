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

    # Check for glucose moiety (simplified hexose ring with 5 hydroxyl groups)
    glucose_pattern = Chem.MolFromSmarts("C1(CO)OC(O)C(O)C(O)C1O")
    glucose_matches = mol.GetSubstructMatches(glucose_pattern)
    if not glucose_matches:
        return False, "No glucose moiety found"

    # Check for ceramide backbone (amide bond connected to long-chain sphingoid base)
    # Amide bond: [NX3][CX3](=O)
    amide_pattern = Chem.MolFromSmarts("N[C](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found"

    # Check for sphingoid base (long-chain amino alcohol)
    # Look for nitrogen atom connected to carbon chain with at least one hydroxyl group
    has_sphingoid = False
    for amide_match in amide_matches:
        amide_nitrogen_idx = amide_match[0]
        amide_nitrogen = mol.GetAtomWithIdx(amide_nitrogen_idx)
        # Get connected atoms to the nitrogen
        neighbors = amide_nitrogen.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Traverse the chain to find an oxygen (hydroxyl group)
                path = Chem.DepthFirstSearch(Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("C(=O)N")))
                has_hydroxyl = False
                atom_idx = neighbor.GetIdx()
                visited = set()
                stack = [atom_idx]
                while stack:
                    idx = stack.pop()
                    if idx in visited:
                        continue
                    visited.add(idx)
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() == 8:  # Oxygen
                        if mol.GetBondBetweenAtoms(atom.GetIdx(), idx).GetBondType() == Chem.rdchem.BondType.SINGLE:
                            has_hydroxyl = True
                            break
                    elif atom.GetAtomicNum() == 6:  # Carbon
                        for nbr in atom.GetNeighbors():
                            if nbr.GetIdx() not in visited:
                                stack.append(nbr.GetIdx())
                if has_hydroxyl:
                    has_sphingoid = True
                    break
        if has_sphingoid:
            break
    if not has_sphingoid:
        return False, "No sphingoid base found"

    # Check for glycosidic linkage between glucose and sphingoid base
    # Oxygen atom connecting glucose ring to sphingoid base
    glycosidic_link_found = False
    for glucose_match in glucose_matches:
        glucose_ring_atoms = set(glucose_match)
        for atom_idx in glucose_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Look for oxygen atoms connected to the glucose ring that are not part of the ring
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if nbr_idx not in glucose_ring_atoms and neighbor.GetAtomicNum() == 8:  # Oxygen
                    # Check if this oxygen connects to the sphingoid base
                    for nbr2 in neighbor.GetNeighbors():
                        if nbr2.GetIdx() != atom_idx and nbr2.GetAtomicNum() == 6:
                            # Check if the connected carbon is part of sphingoid base
                            connected_to_nitrogen = False
                            for path in Chem.FindAllPathsOfLengthN(mol, nbr2.GetIdx(), 4, useBonds=False):
                                if amide_nitrogen_idx in path:
                                    connected_to_nitrogen = True
                                    break
                            if connected_to_nitrogen:
                                glycosidic_link_found = True
                                break
                if glycosidic_link_found:
                    break
            if glycosidic_link_found:
                break
        if glycosidic_link_found:
            break
    if not glycosidic_link_found:
        return False, "No glycosidic linkage between glucose and sphingoid base found"

    # Optional: Check for long hydrocarbon chains (typical of ceramides)
    # Count carbons in acyl chain (fatty acid)
    fatty_acid_length = 0
    for amide_match in amide_matches:
        carbonyl_c_idx = amide_match[1]
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        # Get the carbon chain connected to the carbonyl carbon
        chain = Chem.GetSymmSSSR(mol)
        fatty_acid_carbons = set()
        stack = [neighbor.GetIdx() for neighbor in carbonyl_c.GetNeighbors() if neighbor.GetAtomicNum() == 6]
        while stack:
            idx = stack.pop()
            if idx in fatty_acid_carbons or idx == carbonyl_c_idx:
                continue
            fatty_acid_carbons.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in fatty_acid_carbons:
                        stack.append(nbr.GetIdx())
        fatty_acid_length = max(fatty_acid_length, len(fatty_acid_carbons))
    if fatty_acid_length < 12:
        return False, f"Fatty acid chain too short ({fatty_acid_length} carbons), expected at least 12 carbons"

    return True, "Contains ceramide backbone with glucose attached via glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37547',
                              'name': 'glucosylceramide',
                              'definition': 'Any of the cerebrosides in which the monosaccharide head group is glucose.',
                              'parents': ['CHEBI:33563']},
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
        'success': True,
        'best': True,
        'error': '',
        'stdout': None
    }