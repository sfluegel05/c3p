"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group
    to an aldehydic or ketonic carbonyl group (1-oxacycloalkan-2-ol or an unsaturated analogue).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a lactol, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a comprehensive SMARTS pattern for lactol structure
    # This pattern attempts to capture the hemiacetal structure with oxygen in a ring (1-oxacycle),
    # connected to a hydroxyl group and a carbon which was originally carbonyl (C=O)
    lactol_pattern = Chem.MolFromSmarts('C1CO[C@@H]([OH1])[C@,C](=O)O1')
    
    if mol.HasSubstructMatch(lactol_pattern):
        return True, "Contains a lactol structural pattern (cyclic ether with hydroxyl group and adjacent carbonyl-derived center)"

    return False, "No lactol pattern found (cyclic hemiacetal)"