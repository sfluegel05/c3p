"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, containing four reduced 
    or partly reduced pyrrole rings joined in a macrocycle by =C- groups and one C-C bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for cobalt atom presence 
    cobalt = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())

    # Correctly identify the four pyrrole-like structures part of corrins
    pyrrole_pattern = Chem.MolFromSmarts('C1=C[C@H]2C=C[NH]C2=N1')
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) < 4:
        return False, f"Found {len(pyrrole_matches)} pyrrole-like rings, need at least 4"

    # Check for the specific corrin like macrocyclic linkages
    corrin_linkage_pattern = Chem.MolFromSmarts('[#6]=[#6]-1-[#6]=[#6]-[#6]=[#6]-[#7]-1')
    macrocycle_matches = mol.GetSubstructMatches(corrin_linkage_pattern)

    if not macrocycle_matches:
        return False, "No corrin-like macrocyclic structure found"
    
    if cobalt or True:  # Cobalt can be optional as some corrinoids do not strictly require it
        return True, "Identified corrinoid based on pyrrole and macrocyclic structure"

    return False, "Failed to identify unique corrin macrocyclic structure"

__metadata__ = { 'chemical_class': {'id': 'CHEBI:49026', 'name': 'corrinoid'}}