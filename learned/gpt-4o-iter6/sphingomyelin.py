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
    
    # Define the amide linkage pattern (N-C(=O)-C)
    amide_pattern = Chem.MolFromSmarts("NC(=O)C")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing characteristic amide linkage with fatty acid"
    
    # Define phosphorylcholine group pattern (-COP(=O)(OCC[N+](C)(C)C)-)
    phosphorylcholine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphorylcholine_pattern):
        return False, "Missing or incorrectly structured phosphorylcholine group"
    
    # Define the sphingoid base pattern including necessary stereochemistry and hydroxyl positions
    sphingoid_base_pattern = Chem.MolFromSmarts("[C@H](O)C[C@H](O)")
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "Missing critical hydroxylated sphingoid base features"
    
    return True, "Molecule matches the core structural features of sphingomyelin"