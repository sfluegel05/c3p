"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: CHEBI:XXXXX prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for quinone core (1,4-benzoquinone or 1,4-naphthoquinone)
    quinone_pattern = Chem.MolFromSmarts("[#6]1([#6]=[#8])[#6]=[#6][#6](=[#8])[#6]=[#6]1 |(1:2,3:4,5:6)|")
    if not mol.HasSubstructMatch(quinone_pattern):
        quinone_pattern = Chem.MolFromSmarts("[#6]1([#6]=[#8])[#6]=[#6][#6](=[#8])[#6]=[#6][#6]=[#6]1 |(1:2,3:4,5:6,7:8)|")
        if not mol.HasSubstructMatch(quinone_pattern):
            return False, "No quinone core found"

    # Look for polyprenyl side-chain (long isoprenoid chain)
    # Pattern for isoprene units: C=C(C)C or C=C(C)CC
    isoprene_pattern = Chem.MolFromSmarts("[#6]=[#6](-[#6])-[#6]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, need at least 2"

    # Check if isoprene units are connected in a chain
    # Ensure that the isoprene units are connected sequentially
    isoprene_chain = False
    for match in isoprene_matches:
        for atom in match:
            neighbors = mol.GetAtomWithIdx(atom).GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetIdx() in [a for m in isoprene_matches for a in m]:
                    isoprene_chain = True
                    break
            if isoprene_chain:
                break
        if isoprene_chain:
            break

    if not isoprene_chain:
        return False, "Isoprene units not connected in a chain"

    # Check molecular weight - prenylquinones typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for prenylquinone"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for prenylquinone"
    if o_count < 2:
        return False, "Must have at least 2 oxygens (quinone core)"

    return True, "Contains quinone core with polyprenyl side-chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:XXXXX',
                          'name': 'prenylquinone',
                          'definition': 'A quinone substituted by a polyprenyl-derived side-chain. Prenylquinones occur in all living cells. Due to their amphiphilic character, they are mainly located in biological membranes where they function as electron and proton carriers in the photosynthetic and respiratory electron transport chains. Some prenylquinones also perform more specialised roles sucy as antioxidants and enzyme cofactors. Prenylquinones are classified according to ring structure: the main classes are menaquinones, phylloquinones, ubiquinones and plastoquinones.'},
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
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}