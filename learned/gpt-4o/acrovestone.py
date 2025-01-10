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

    # Broadly define isoflavonoid backbone pattern
    isoflavone_pattern = Chem.MolFromSmarts("c1cc2oc(ccc2c1)-c3ccc(O)c(=O)c3")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavonoid polyphenolic core detected"
    
    # Define diverse glycoside patterns
    glycoside_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1"), # Glucosides
        Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"),         # Fucosides, rhamnosides
        Chem.MolFromSmarts("O1C([C@@H](O)[C@H](O)C(O)C1)CO"),                 # Misc hexoses
    ]
    
    # Check for glycosidic attachments
    if not any(mol.HasSubstructMatch(pattern) for pattern in glycoside_patterns):
        return False, "No glycoside moieties detected"

    # Verify presence of polyphenolic structure: at least three aromatic rings
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 3:
        return False, "Too few aromatic rings for a polyphenolic structure"

    # Check for characteristic methoxy (-OCH3) or hydroxyl (-OH) substitutions
    methoxy_pattern = Chem.MolFromSmarts("COC")
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    
    if not (mol.HasSubstructMatch(methoxy_pattern) or mol.HasSubstructMatch(hydroxyl_pattern)):
        return False, "Methoxy or hydroxyl group substitutions missing from aromatic rings"

    # All checks passed
    return True, "Molecule is consistent with structural features of acrovestone"