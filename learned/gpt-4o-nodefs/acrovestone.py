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
        Chem.MolFromSmarts("O=C1C2=C(OC=C1)C=CC=C2"),  # Standard chromen structure
        Chem.MolFromSmarts("C1=C2C(=O)C=CC=C2C=CC2=C1C=CC=C2"),  # Extended aromatic forms
        Chem.MolFromSmarts("O=C1C=CC2=CC=C(O)C=C2O1"),  # More generic isoflavone structures
        # Add more extensive patterns if necessary
    ]

    # Ensure the presence of one isoflavone core structure 
    if not any(mol.HasSubstructMatch(pat) for pat in isoflavone_patterns):
        return False, "No isoflavone core structure found"
    
    # Check for glycosidic patterns (Generic pattern including potential O-glycosidic bonds)
    glycoside_patterns = [
        Chem.MolFromSmarts("COC1C(O)C(O)C(O)C(O)C1O"),  # Common glucose ring
        Chem.MolFromSmarts("COC1C(C(O)C(C(C1O)O)O)O"),  # Rhamnose or other sugar variants
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"),  # Generic glycoside with stereochemistry
        # Additional saccharide attachments
    ]
    
    # Verify a matching glycosidic pattern
    if not any(mol.HasSubstructMatch(pat) for pat in glycoside_patterns):
        return False, "No recognizable glycosidic bonds (sugar moieties) attached"
    
    return True, "Typical acrovestone-like structure detected: isoflavone with glycosidic linkage"