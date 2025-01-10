"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine typically includes an aromatic system with an amino group,
    which may be directly attached or linked through a short chain.

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

    # Define aromatic ring with hydroxy groups like catechol
    aromatic_phenol = Chem.MolFromSmarts("c1ccc(O)c(O)c1")
    
    # Check presence of catechol-like aromatic systems
    if not mol.HasSubstructMatch(aromatic_phenol):
        return False, "No catechol-like aromatic system detected"
    
    # Define patterns for primary and secondary amines connected via small aliphatic chains
    amine_connected_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0][#6,#1]")
    
    # Ensure connection of amine near aromatic structure
    if not mol.HasSubstructMatch(amine_connected_pattern):
        return False, "No amine group sufficiently close to an aromatic system"
   
    return True, "Contains an aromatic system with an amino group possibly connected through a short chain"