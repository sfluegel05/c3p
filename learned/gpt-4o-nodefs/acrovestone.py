"""
Classifies: CHEBI:2440 acrovestone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone-like structure based on its SMILES string.
    An acrovestone-like molecule typically contains an isoflavone core bonded to a glycosyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches acrovestone-like structures, False otherwise
        str: Reason for classification or lack thereof
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expanded pattern for isoflavone core (4H-chromen-4-one derivatives)
    isoflavone_patterns = [
        Chem.MolFromSmarts("O=C1C=C(O)C2=CC=CC=C2C1C"),  # main flavone pattern
        Chem.MolFromSmarts("O=C1C=C(O)C=CC=C2C1=CC=CC2") # more generic isoflavone structure
    ]

    # Ensure the presence of one isoflavone core structure 
    if not any(mol.HasSubstructMatch(pat) for pat in isoflavone_patterns):
        return False, "No isoflavone core structure found"
    
    # Check for glycosidic patterns (Generic pattern including potential O-glycosidic  bonds)
    glycoside_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@@H]([C@H](O)[C@@H](O)[C@H]1O)"),  # Polyol Ring
        Chem.MolFromSmarts("O[C@H]1[C@@H]([C@H](O)[C@@H](O)[C@H](O)[C@H]1O)")  # Variations of glycosides
        # Add patterns for common sugar attachments- both linear and cyclic forms
    ]
    
    # Verify a matching glycosidic pattern
    if not any(mol.HasSubstructMatch(pat) for pat in glycoside_patterns):
        return False, "No recognizable glycosidic bonds (sugar moieties) attached"
    
    return True, "Typical acrovestone-like structure detected: isoflavone with glycosidic linkage"