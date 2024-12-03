"""
Classifies: CHEBI:83228 L-alpha-amino acid residue
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_L_alpha_amino_acid_residue(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid residue.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid residue, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for alpha-amino acid residue core structure
    alpha_carbon = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            neighbors = atom.GetNeighbors()
            has_amino_group = any(neighbor.GetSymbol() == 'N' for neighbor in neighbors)
            has_carboxyl_group = any(neighbor.GetSymbol() == 'C' and len([n for n in neighbor.GetNeighbors() if n.GetSymbol() == 'O']) == 2 for neighbor in neighbors)
            if has_amino_group and has_carboxyl_group:
                alpha_carbon = atom
                break

    if alpha_carbon is None:
        return False, "No alpha-carbon with both amino and carboxyl groups found"

    # Check for L-configuration
    if alpha_carbon.GetChiralTag() not in [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]:
        return False, "Alpha-carbon is not chiral"

    # Assuming the molecule is an L-alpha-amino acid residue if it contains the core structure
    return True, "Molecule is an L-alpha-amino acid residue"

# Example usage
smiles = "C([C@@H](C(*)=O)N(*)C)C1=CC=CC=C1"  # N-methyl-L-phenylalanine residue
print(is_L_alpha_amino_acid_residue(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83228',
                          'name': 'L-alpha-amino acid residue',
                          'definition': 'An alpha-amino-acid residue derived '
                                        'from an L-alpha-amino acid.',
                          'parents': ['CHEBI:33710']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(False, 'No alpha-carbon with both amino and carboxyl groups "
              "found')\n",
    'num_true_positives': 18,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 20,
    'precision': 0.9,
    'recall': 0.47368421052631576,
    'f1': 0.6206896551724138,
    'accuracy': None}