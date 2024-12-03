"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of both amino and carboxyl groups
    amino_group = False
    carboxyl_group = False
    alpha_carbon = None

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetTotalNumHs() == 2:
            amino_group = True
            alpha_carbon = atom.GetNeighbors()[0]  # alpha carbon is the carbon bonded to the amino group
        if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
                    carboxyl_group = True

    if not amino_group or not carboxyl_group:
        return False, "Missing amino or carboxyl group"

    # Check for L-configuration at the alpha-carbon
    if alpha_carbon is None:
        return False, "Alpha carbon not found"

    chiral_tag = alpha_carbon.GetChiralTag()
    if chiral_tag not in [Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW]:
        return False, "Alpha carbon is not chiral"

    # Ensure the alpha carbon is bonded to the carboxyl group and another carbon
    alpha_carbon_bonds = [neighbor.GetSymbol() for neighbor in alpha_carbon.GetNeighbors()]
    if alpha_carbon_bonds.count('C') < 2 or alpha_carbon_bonds.count('O') < 1:
        return False, "Alpha carbon does not have the correct bonding pattern"

    # Check if the configuration is L
    if alpha_carbon.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
        return True, "L-alpha-amino acid with CW configuration"
    elif alpha_carbon.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return True, "L-alpha-amino acid with CCW configuration"
    else:
        return False, "Alpha carbon has an unknown chiral configuration"

    return False, "Does not meet L-alpha-amino acid criteria"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15705',
                          'name': 'L-alpha-amino acid',
                          'definition': 'Any alpha-amino acid having '
                                        'L-configuration at the alpha-carbon.',
                          'parents': ['CHEBI:33704']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 44,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}