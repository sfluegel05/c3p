"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: CHEBI:25657 3'-hydroxyflavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Tuple

def is_3__hydroxyflavanones(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone is a flavanone with a hydroxy substituent at position 3' of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enumerate tautomers and stereoisomers
    tautomers = [Chem.MolFromSmiles(Chem.MolToSmiles(taut, isomericSmiles=True)) for taut in AllChem.EnumerateTautomers(mol)]
    stereoisomers = list({Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(isomer, isomericSmiles=True)), isomericSmiles=True)
                          for isomer in AllChem.EnumerateStereoisomers(mol)})
    
    # Check if any tautomer/stereoisomer matches the 3'-hydroxyflavanone pattern
    for tautomer in tautomers:
        for stereoisomer in stereoisomers:
            if _is_3__hydroxyflavanone(stereoisomer):
                return True, "Matches the 3'-hydroxyflavanone structural pattern"
    
    return False, "Does not match the 3'-hydroxyflavanone structural pattern"

def _is_3__hydroxyflavanone(mol: Chem.Mol) -> bool:
    """
    Checks if a molecule matches the 3'-hydroxyflavanone structural pattern.

    Args:
        mol (Chem.Mol): RDKit molecule object

    Returns:
        bool: True if the molecule matches the 3'-hydroxyflavanone pattern, False otherwise
    """
    
    # Check for flavanone scaffold
    flavanone_pattern = Chem.MolFromSmarts("[O;R]c1ccc(cc1)[C@H]2CC(=O)Oc2")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False
    
    # Check for hydroxy group at 3' position
    hydroxy_3__pattern = Chem.MolFromSmarts("c1ccc(cc1O)c2ccc(cc2)")
    if not mol.HasSubstructMatch(hydroxy_3__pattern):
        return False
    
    # Additional checks for common substituents (prenyl, glycosidic)
    prenyl_pattern = Chem.MolFromSmarts("[CH2]=[CH][CH2]=[CH][CH3]")
    glycosidic_pattern = Chem.MolFromSmarts("[OX2]C[C@@H]([C@H]([C@@H]([C@H](O[C@@H]1[C@H]([C@@H]([C@@H](O1)O)O)O)O)O)O")
    
    if mol.HasSubstructMatch(prenyl_pattern) or mol.HasSubstructMatch(glycosidic_pattern):
        return True
    
    return True