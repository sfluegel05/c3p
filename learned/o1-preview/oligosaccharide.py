"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is composed of 2 to 20 monosaccharide units joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize variables
    sugar_rings = []
    glycosidic_bonds = []

    # Detect rings in the molecule
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    # Identify sugar rings (5 or 6 membered rings with one oxygen)
    for ring in atom_rings:
        if len(ring) in [5,6]:
            o_count = 0
            c_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    o_count += 1
                elif atom.GetAtomicNum() == 6:
                    c_count += 1
            # Sugar rings have exactly one oxygen in the ring
            if o_count == 1 and c_count == (len(ring)-1):
                sugar_rings.append(set(ring))

    num_monosaccharides = len(sugar_rings)

    if num_monosaccharides < 2:
        return False, f"Found {num_monosaccharides} monosaccharide unit(s), need at least 2 for oligosaccharide"

    # Identify glycosidic linkages
    for bond in mol.GetBonds():
        # Look for bonds connecting two sugar rings via an oxygen
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 8 or a2.GetAtomicNum() == 8:
            # One atom is oxygen
            other_atom = a1 if a2.GetAtomicNum() == 8 else a2
            o_atom = a2 if a2.GetAtomicNum() == 8 else a1
            # Oxygen connected to carbons in two different sugar rings
            if other_atom.GetAtomicNum() == 6:
                # Check if the carbon atom is part of a sugar ring
                for ring in sugar_rings:
                    if other_atom.GetIdx() in ring:
                        # Check if oxygen is connected to another carbon in a different ring
                        for nbr in o_atom.GetNeighbors():
                            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != other_atom.GetIdx():
                                for other_ring in sugar_rings:
                                    if nbr.GetIdx() in other_ring and other_ring != ring:
                                        glycosidic_bonds.append(bond)
    num_glycosidic_bonds = len(glycosidic_bonds)

    if num_glycosidic_bonds < 1:
        return False, "No glycosidic linkages found between monosaccharide units"

    # Check if the number of glycosidic bonds corresponds to number of monosaccharides minus one
    if num_glycosidic_bonds < num_monosaccharides - 1:
        return False, f"Insufficient glycosidic linkages: found {num_glycosidic_bonds}, expected at least {num_monosaccharides - 1}"

    # Additional check: ensure that the molecule is not too large (exclude polysaccharides)
    if num_monosaccharides > 20:
        return False, f"Found {num_monosaccharides} monosaccharide units, which may indicate a polysaccharide"

    return True, f"Contains {num_monosaccharides} monosaccharide units connected via glycosidic linkages"

__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'oligosaccharide',
                              'definition': 'A compound in which monosaccharide units are joined by glycosidic linkages. The term is commonly used to refer to a defined structure as opposed to a polymer of unspecified length or a homologous mixture. When the linkages are of other types the compounds are regarded as oligosaccharide analogues.',
                              'parents': []},
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
        'success': None,
        'best': None,
        'error': '',
        'stdout': None,
        'num_true_positives': None,
        'num_false_positives': None,
        'num_true_negatives': None,
        'num_false_negatives': None,
        'num_negatives': None,
        'precision': None,
        'recall': None,
        'f1': None,
        'accuracy': None}