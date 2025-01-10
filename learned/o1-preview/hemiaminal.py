"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is any organic amino compound that has an amino group and a hydroxy group attached to the same carbon atom.
    Hemiaminals are intermediates in the formation of imines by addition of an amine to an aldehyde or ketone.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a carbon with both hydroxyl and amino groups attached,
    # excluding carbonyls (C=O), imines (C=N), and nitriles (C#N)
    hemiaminal_pattern = Chem.MolFromSmarts('[C;!$([C]=[O]);!$([C]=[N]);!$([C]#[N])](O)(N)')
    if hemiaminal_pattern is None:
        return None, "Invalid SMARTS pattern for hemiaminal"
    
    # Check for hemiaminal substructure
    if mol.HasSubstructMatch(hemiaminal_pattern):
        return True, "Molecule contains a hemiaminal group (carbon attached to both -OH and -NH2/-NHR/-NR2)"
    else:
        return False, "No hemiaminal group found"

__metadata__ = {
    'chemical_class': {   
        'name': 'hemiaminal',
        'definition': 'Any organic amino compound that has an amino group and a hydroxy group attached to the same carbon atom. Hemiaminals are intermediates in the formation of imines by addition of an amine to an aldehyde or ketone; those derived from primary amines are particularly unstable.',
    }
}