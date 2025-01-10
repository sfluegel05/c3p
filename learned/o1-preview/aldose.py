"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: aldose
"""

from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is an aldehydic sugar (polyhydroxy aldehyde) or its intramolecular hemiacetal form
    (a cyclic structure like a furanose or pyranose ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Too few carbons for an aldose"

    # Check for aldehyde group (open-chain form)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H](=O)")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)

    # Count the number of hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    num_hydroxyl = len(mol.GetSubstructMatches(hydroxyl_pattern))

    if has_aldehyde:
        if num_hydroxyl >= 2:
            return True, f"Contains aldehyde group and {num_hydroxyl} hydroxyl groups"
        else:
            return False, "Contains aldehyde group but insufficient hydroxyl groups"
    else:
        # Check for cyclic hemiacetal forms (furanose or pyranose rings)
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        found_hemiacetal = False

        for ring in atom_rings:
            ring_size = len(ring)
            if ring_size == 5 or ring_size == 6:
                num_oxygen = 0
                num_carbons = 0
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() == 8:
                        num_oxygen += 1
                    elif atom.GetAtomicNum() == 6:
                        num_carbons += 1
                if ((ring_size == 5 and num_carbons == 4 and num_oxygen == 1) or
                    (ring_size == 6 and num_carbons == 5 and num_oxygen == 1)):
                    found_hemiacetal = True
                    break

        if found_hemiacetal:
            if num_hydroxyl >= 2:
                return True, f"Contains cyclic hemiacetal (furanose or pyranose ring) and {num_hydroxyl} hydroxyl groups"
            else:
                return False, "Contains cyclic hemiacetal but insufficient hydroxyl groups"
        else:
            return False, "Does not contain aldehyde group or cyclic hemiacetal"

__metadata__ = {   'chemical_class': {   'name': 'aldose',
                              'definition': 'Aldehydic parent sugars (polyhydroxy '
                                            'aldehydes H[CH(OH)]nC(=O)H, n >= 2) '
                                            'and their intramolecular hemiacetals.'},
        'config': {   'llm_model_name': 'your_llm_model_name_here',
                      'f1_threshold': 0.8,
                      'max_attempts': 5},
        'message': None,
        'success': True}