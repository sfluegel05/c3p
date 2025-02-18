"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids contain an aromatic ring connected to a propane chain or similar structures
    like flavonoids, coumarins, etc.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expand aromatic systems 
    aromatic_system_patterns = [
        Chem.MolFromSmarts("c1ccccc1"),  # Benzene
        Chem.MolFromSmarts("c1ccc2c(c1)cccc2"),  # Naphthalene
        Chem.MolFromSmarts("c1cc(c2cocc2)ccc1"),  # Coumarin
        Chem.MolFromSmarts("c1c2c(ccc1)c(=O)oc2"),  # Flavonoids basic core
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in aromatic_system_patterns):
        return False, "No suitable aromatic system found"

    # Improved phenylpropanoid backbone (simplified, generic)
    phenylpropanoid_backbone_patterns = [
        Chem.MolFromSmarts("c-[CX3](=O)-O"),  # Coumarin-related ester linkage
        Chem.MolFromSmarts("c1cc(c2c(c1)ccc(c2)O)"),  # Additional flavonoid-like
        Chem.MolFromSmarts("c-C-C(=O)-OC"),  # Basic ester backbone
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in phenylpropanoid_backbone_patterns):
        return False, "No recognizable phenylpropanoid backbone found"

    # More comprehensive functional groups covering phenylpropanoids
    functional_group_patterns = [
        Chem.MolFromSmarts("c-[OH]"),      # Phenolic OH groups
        Chem.MolFromSmarts("c-CO"),        # Carbonyl
        Chem.MolFromSmarts("c-COC"),       # Methoxy or similar groups
        Chem.MolFromSmarts("c=C[OX2]"),    # Enol or ester groups
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in functional_group_patterns):
        return False, "No common phenylpropanoid functional groups found"

    return True, "Structure contains an aromatic system and functional groups characteristic of a phenylpropanoid."