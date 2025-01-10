"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is an aromatic amide obtained by acylation of aniline,
    which typically includes a phenyl group directly attached to a nitrogen of an amide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Anilide pattern: phenyl group directly bonded to amide nitrogen
    # We will use a pattern that allows for general anilide structure
    # Phenyl ring (c1ccccc1) connected to amide nitrogen (NC(=O))
    # This pattern is broad to capture variations, focuses on the critical aspect of the anilide structure.
    
    anilide_pattern = Chem.MolFromSmarts("c1ccccc1NC(=O)")
    
    # Check for anilide pattern with flexibility in C part of NC(=O)
    match = mol.HasSubstructMatch(anilide_pattern)
    
    if not match:
        return False, "No anilide pattern found (phenyl group directly bonded to amide nitrogen)"

    return True, "Contains phenyl group directly bonded to amide nitrogen, fitting anilide definition"