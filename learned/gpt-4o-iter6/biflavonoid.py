"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem

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
    # Parse SMILES using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refine SMARTS pattern for flavonoid-like structures
    # More generic flavonoid unit pattern to account for diversity:
    pyran_pattern = Chem.MolFromSmarts("c1c(O)cc(cc1)-C2=O")  # A basic Acceptor-Matching pattern in flavonoids
    
    # Finding matches for each flavonoid-like substructure
    flavonoid_like_units = mol.GetSubstructMatches(pyran_pattern)

    # Check for presence of at least two such units
    if len(flavonoid_like_units) < 2:
        return False, "Less than two flavonoid units identified"
    
    # Check for potential coupling linkages - a single atom/bond connection:
    # Loop through combinations
    connection_found = False
    for i in range(len(flavonoid_like_units)):
        for j in range(i+1, len(flavonoid_like_units)):
            # check adjacency and connection oxygen/carbon
            adj_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in ["C", "O"]]
            if any(
                mol.GetBondBetweenAtoms(adj_atoms[k].GetIdx(), adj_atoms[i].GetIdx()) 
                for k in range(len(adj_atoms))
            ):
                connection_found = True
                break

    if not connection_found:
        return False, "No connecting single atom or bond found between flavonoid units"
    
    return True, "Contains at least two flavonoid units linked by a single bond or atom"

# Metadata remains the same as the initial attempt
__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:biflavonoid',
        'name': 'biflavonoid',
        'definition': 'A flavonoid oligomer that is obtained by the oxidative coupling of at least two units of aryl-substituted benzopyran rings or its substituted derivatives, resulting in the two ring systems being joined together by a single atom or bond.',
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