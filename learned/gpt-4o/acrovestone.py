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

    # 1. Check for polyphenolic structure: look for isoflavonoid backbones
    # Isoflavones: characteristic C3-C2-C1=O structure with substituted aromatic rings
    isoflavone_pattern = Chem.MolFromSmarts("c1cc(O)c(C(=O)c2coc(c12)c3ccc(OC)c(O)c3)")
    
    # Check for any isoflavonoid core
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavonoid polyphenolic core detected"

    # 2. Expanded glycosylation patterns: look for various sugar moieties
    glycoside_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"),  # Glucosides
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO1"),     # Other hexoses
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)CO1")              # Pentoses
    ]
    
    # Check for any glycosidic attachments
    if not any(mol.HasSubstructMatch(pattern) for pattern in glycoside_patterns):
        return False, "No glycoside moieties detected"

    # 3. Verify presence of adequate number of aromatic rings which could imply an extended polyphenol
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 3:
        return False, "Too few aromatic rings for a polyphenolic structure"

    # All checks passed
    return True, "Molecule is consistent with structural features of acrovestone"