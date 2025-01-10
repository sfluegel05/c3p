"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is a N-hydroxy-alpha-amino-acid based on its SMILES string.
    This requires an amino acid in which at least one hydrogen attached to the amino group 
    is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for an alpha amino acid backbone
    # This is more generic, allowing any substitution and focusing on the carboxy and amino group pattern
    amino_acid_pattern = Chem.MolFromSmarts("C(C(=O)O)[NH2,NH,N]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No alpha amino acid backbone found"

    # Detect if there is an N-Hydroxy group attached to any nitrogen in the substructure
    has_n_hydroxy = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check for nitrogen atoms
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:  # Check for an attached oxygen
                    has_n_hydroxy = True
                    break
    
    if not has_n_hydroxy:
        return False, "No N-hydroxy modification found"

    return True, "Contains amino acid backbone with N-hydroxy modification"

__metadata__ = {  
    'chemical_class': {   
        'name': 'N-hydroxy-alpha-amino-acid',
        'definition': 'Any amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.',
    },
    'message': None,
    'attempt': 0,
    'success': True,
}