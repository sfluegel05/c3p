"""
Classifies: CHEBI:2440 acrovestone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is of the acrovestone class.
    Acrovestone is characterized as a polyphenolic compound, typically involving isoflavonoid cores 
    and various glycoside attachments, often with methoxy and hydroxyl substitutions on aromatic rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as acrovestone, False otherwise
        str: Reason for classification
    """
    
    # Try to parse the SMILES string, check validity
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Revised definition for isoflavonoid backbone pattern
    isoflavone_pattern = Chem.MolFromSmarts("c1c(ccc2c1oc(=O)c3ccccc3-2)O")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavonoid polyphenolic core detected"
    
    # Expanded glycoside patterns to capture a wider array of possible configurations
    glycoside_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1[*]"), # Various hexoses
        Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"),         # Fucosides, rhamnosides
        Chem.MolFromSmarts("O1[C@@H]([C@H](O)[C@H](O)C1)CO"),                 # Misc hexoses
    ]
    
    # Check for glycosidic attachments
    if not any(mol.HasSubstructMatch(pattern) for pattern in glycoside_patterns):
        return False, "No glycoside moieties detected"

    # Verify presence of at least two aromatic rings for polyphenolic structure
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 2:
        return False, "Too few aromatic rings for a polyphenolic structure"

    # Check for characteristic methoxy or hydroxyl substitutions on aromatic rings
    methoxy_pattern = Chem.MolFromSmarts("c-COC")
    hydroxyl_pattern = Chem.MolFromSmarts("c-O")
    
    if not (mol.HasSubstructMatch(methoxy_pattern) or mol.HasSubstructMatch(hydroxyl_pattern)):
        return False, "Methoxy or hydroxyl group substitutions missing from aromatic rings"

    # All checks passed
    return True, "Molecule is consistent with structural features of acrovestone"