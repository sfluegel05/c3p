"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: Alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    Alkaloids are typically nitrogen-containing organic compounds with complex and diverse structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely an alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of nitrogen atoms
    num_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if num_nitrogens == 0:
        return False, "No nitrogen atoms found"

    # Check for ring structures
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized() or ring_info.NumRings() < 1:
        return False, "No ring structures found"
        
    # Additional checks could be added here based on common alkaloid patterns

    return True, "Contains nitrogen atoms and ring structures, indicating possible alkaloid structure"

__metadata__ = {   'chemical_class': {   'name': 'alkaloid',
                          'definition': 'Alkaloids are naturally occurring '
                                        'organic compounds that mostly contain '
                                        'basic nitrogen atoms.'}
}