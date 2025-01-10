"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid contains an indole skeleton and typically includes nitrogen within a polycyclic framework.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexible pattern to identify an indole or indole-derived structure
    # Without relying strictly on the full aromatic ring.
    indole_patterns = [
        Chem.MolFromSmarts("n1c2ccccc2c1"),  # general indole, ignoring aromaticity
        Chem.MolFromSmarts("c1cnc2ccccc2c1"),  # another variant
    ]
    
    # Check against each pattern
    if not any(mol.HasSubstructMatch(pattern) for pattern in indole_patterns):
        return False, "No recognizable indole-like skeleton found"

    # Check for presence of nitrogen beyond the basic indole nitrogen
    # Alkaloids may have even integrated nitrogen in non-basic ring forms
    extra_nitrogen_count = sum(1 for atom in mol.GetAtoms() 
                               if atom.GetAtomicNum() == 7 and atom.GetDegree() > 1)
    
    if extra_nitrogen_count < 1:
        return False, "No extra nitrogen atoms indicating possible alkaloid structure"

    return True, "Contains a flexible indole-like skeleton and potential alkaloid features"