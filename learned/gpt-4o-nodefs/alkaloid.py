"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: Alkaloid
"""
from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    Alkaloids are typically nitrogen-containing organic compounds with complex and diverse structures,
    often containing nitrogen atoms within heterocyclic rings.

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

    # Check for nitrogen atom in a ring
    has_nitrogen_in_ring = False
    for nitrogen in [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]:
        if nitrogen.IsInRing():
            has_nitrogen_in_ring = True
            break

    if not has_nitrogen_in_ring:
        return False, "Nitrogen not found in ring structures"

    # Additional pattern checks can be done here, such as for common alkaloid ring structures

    return True, "Contains nitrogen atoms in ring structures, indicating possible alkaloid structure"

__metadata__ = {   'chemical_class': {   'name': 'alkaloid',
                          'definition': 'Alkaloids are naturally occurring '
                                        'organic compounds that mostly contain '
                                        'basic nitrogen atoms.'}
}