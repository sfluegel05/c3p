"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine typically contains an amino group linked to an aromatic system,
    sometimes connected directly or through a short spacer.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES 
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General aromatic system pattern
    aromatic_pattern = Chem.MolFromSmarts("c")  # Any aromatic carbon

    # Check for aromatic system
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic system (with carbon) found"
    
    # Look for a generalized amine pattern linked possibly close to the aromatic system; relaxed to allow for zero or more atoms in the middle
    # [NX3;H2,H1,H0] matches any nitrogen (primary, secondary, tertiary amines)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]C*[*]c")  # Relaxing linkage requirements
    
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine involving an aromatic system and short connectors found"

    return True, "Contains an aromatic system with an amino group potentially connected through a shorter chain"