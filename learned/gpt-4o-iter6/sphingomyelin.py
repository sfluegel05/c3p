"""
Classifies: CHEBI:64583 sphingomyelin
"""
from rdkit import Chem

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string and check for valid molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General pattern for sphingoid base (hydroxylated long chain with possible amine)
    sphingoid_base_pattern = Chem.MolFromSmarts("[CX3,CX4][CX2,CX3,CX4](O)[CX2,CX3,CX4][CX2,CX3](N)")
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "Missing or incomplete critical sphingoid base features"

    # Pattern for amide linkage, expanded to capture context as necessary
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing characteristic amide linkage with fatty acid"
    
    # Phosphorylcholine group pattern, allowing more structural flexibility 
    phosphorylcholine_pattern = Chem.MolFromSmarts("C[O,P]([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphorylcholine_pattern):
        return False, "Missing or incorrectly structured phosphorylcholine group"
    
    return True, "Molecule matches the core structural features of sphingomyelin"