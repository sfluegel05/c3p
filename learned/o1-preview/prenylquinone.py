"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

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
    
    # Kekulize molecules to standardize representations
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except:
        pass  # Kekulization can fail for some molecules
    
    # Define a general quinone pattern: ring with two conjugated ketone groups
    quinone_pattern = Chem.MolFromSmarts('[#6]=O[#6]:[#6]:[#6]:[#6]=O')  # General quinone pattern
    
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "No quinone ring found"
    
    # Get the indices of the atoms matching the quinone pattern
    quinone_matches = mol.GetSubstructMatches(quinone_pattern)
    if not quinone_matches:
        return False, "No quinone ring found"
    
    # Collect quinone ring atoms
    quinone_atoms = set()
    for match in quinone_matches:
        quinone_atoms.update(match)
    
    # Define isoprene unit pattern
    isoprene_pattern = Chem.MolFromSmarts('C(=C)C(C)C')  # Isoprene unit
    
    # Identify prenyl side chains attached to quinone ring
    prenyl_side_chain_found = False
    for atom_idx in quinone_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in quinone_atoms:
                # Potential side chain
                side_chain_atoms = set()
                side_chain_bonds = set()
                atoms_to_visit = [neighbor_idx]
                while atoms_to_visit:
                    current_idx = atoms_to_visit.pop()
                    if current_idx not in side_chain_atoms and current_idx not in quinone_atoms:
                        side_chain_atoms.add(current_idx)
                        current_atom = mol.GetAtomWithIdx(current_idx)
                        for nbr in current_atom.GetNeighbors():
                            nbr_idx = nbr.GetIdx()
                            if nbr_idx not in side_chain_atoms and nbr_idx not in quinone_atoms:
                                atoms_to_visit.append(nbr_idx)
                                bond = mol.GetBondBetweenAtoms(current_idx, nbr_idx)
                                side_chain_bonds.add(bond.GetIdx())
                # Create side chain mol
                side_chain = Chem.PathToSubmol(mol, side_chain_bonds)
                # Check for isoprene units
                isoprene_matches = side_chain.GetSubstructMatches(isoprene_pattern)
                if isoprene_matches:
                    prenyl_side_chain_found = True
                    break
        if prenyl_side_chain_found:
            break
    
    if not prenyl_side_chain_found:
        return False, "No prenyl side chain attached to quinone ring"
    
    return True, "Contains quinone ring with polyprenyl-derived side chain"

__metadata__ = {   'chemical_class': {   'id': '',
                              'name': 'prenylquinone',
                              'definition': 'A quinone substituted by a polyprenyl-derived side-chain. Prenylquinones occur in all living cells. Due to their amphiphilic character, they are mainly located in biological membranes where they function as electron and proton carriers in the photosynthetic and respiratory electron transport chains. Some prenylquinones also perform more specialised roles such as antioxidants and enzyme cofactors. Prenylquinones are classified according to ring structure: the main classes are menaquinones, phylloquinones, ubiquinones and plastoquinones.',
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
        'attempt': 1,
        'success': True,
        'best': True,
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