"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is defined as a flavonoid oligomer obtained by the oxidative coupling 
    of at least two units of aryl-substituted benzopyran rings or its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify two benzopyran-like structures within the molecule
    benzopyran_pattern = Chem.MolFromSmarts("c1ccccc1C2=CC=CC=C2O") # Basic benzopyran pattern
    flavonoid_units = mol.GetSubstructMatches(benzopyran_pattern)
    if len(flavonoid_units) < 2:
        return False, "Less than two flavonoid units identified"

    # Check for a single atom or bond linkage between flavonoid units
    # Check for any common atoms or bonds linking flavonoid units
    graph = Chem.GetShortestPathMatrix(mol)
    linking_bonds = 0
    for i in range(len(flavonoid_units)):
        for j in range(i+1, len(flavonoid_units)):
            if graph[flavonoid_units[i][0]][flavonoid_units[j][0]] == 1:
                linking_bonds += 1
    
    if linking_bonds < 1:
        return False, "No single atom or bond link found between flavonoid units"

    return True, "Contains at least two flavonoid units linked by a single bond"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:biflavonoid',
        'name': 'biflavonoid',
        'definition': 'A flavonoid oligomer that is obtained by the oxidative coupling of at least two units of aryl-substituted benzopyran rings or its substituted derivatives, resulting in the two ring systems being joined together by a single atom or bond.',
        'parents': ['CHEBI:'],
    },
    'config': {   
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
}