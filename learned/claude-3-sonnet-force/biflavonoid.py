"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: CHEBI:33136 biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is a flavonoid oligomer obtained by the oxidative coupling of at least two units
    of aryl-substituted benzopyran rings or its substituted derivatives, resulting in the two ring
    systems being joined together by a single atom or bond.

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
    
    # Look for flavonoid subunits
    flavonoid_pattern = Chem.MolFromSmarts("c1c(O)cc2c(c1)C(=O)c3ccccc3O2")
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_pattern)
    
    # Check for at least 2 unique flavonoid subunits
    unique_subunits = set(map(tuple, flavonoid_matches))
    if len(unique_subunits) < 2:
        return False, "Molecule does not contain at least 2 unique flavonoid subunits"
    
    # Check if flavonoid subunits are joined by a single atom or bond
    fragments = Chem.GetAtomJoinedMolecules(mol)
    if len(fragments) == 2:
        fragment1_atoms = set(fragments[0].GetAtoms())
        fragment2_atoms = set(fragments[1].GetAtoms())
        
        # Check if fragments correspond to flavonoid subunits
        if any(set(match).issubset(fragment1_atoms) for match in unique_subunits) and \
           any(set(match).issubset(fragment2_atoms) for match in unique_subunits):
            return True, "Contains at least 2 flavonoid subunits joined by a single atom or bond"
    
    return False, "Flavonoid subunits not joined by a single atom or bond"