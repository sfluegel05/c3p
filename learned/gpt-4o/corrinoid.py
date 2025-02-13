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
    
    # Check for cobalt atom presence (indicator for some corrinoids)
    cobalt = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())

    # SMARTS for reduced/partly reduced pyrrole - adjust to include typical bonds of corrins
    # For simplicity, their pyrrolic nature and typical connectivity is considered
    pyrrole_pattern = Chem.MolFromSmarts('C1=CC=[NH]C=C1')  # Detects pyrrole-like structures
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) < 4:
        return False, f"Found {len(pyrrole_matches)} pyrrole-like rings, need at least 4"
    
    # Additional consideration: further patterns for partial pyrrole reduction, connected via expected corrin linkages
    partial_reduction_pattern = Chem.MolFromSmarts('C1=CC[NH]=C1')  # Simple representation of partial reduction
    reduction_matches = mol.GetSubstructMatches(partial_reduction_pattern)

    # Ensure macrocycle coherence by using a simplified accepted linkage pattern
    # Generally encompasses directly connected =C- groups and a direct C-C characteristic bond present
    macrocycle_confirms = False
    # Typically, manual inspection of node-link relationships (edges) across atoms is preferred here
    macrocycle_pattern = Chem.MolFromSmarts('C(-C=N)-C=N-C=C-C1=N-C=C1')  # Placeholder for typical corrinoid link
    macrocycle_confirmations = mol.GetSubstructMatches(macrocycle_pattern)
    
    if macrocycle_confirmations:
        macrocycle_confirms = True

    if macrocycle_confirms and (cobalt or True):  # Cobalt can be optional
        return True, "Identified corrinoid macrocyclic nucleus with reduced pyrrole rings"
    
    return False, "Failed to identify unique corrin macrocyclic nucleus"

__metadata__ = { 'chemical_class': {'id': 'CHEBI:49026', 'name': 'corrinoid'}}