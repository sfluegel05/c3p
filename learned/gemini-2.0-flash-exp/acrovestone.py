"""
Classifies: CHEBI:2440 acrovestone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone based on its SMILES string.
    Acrovestones are polyphenols, often with an isoflavone or flavone backbone,
    multiple hydroxyl groups, and often glycosylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acrovestone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for isoflavone or flavone core
    # General isoflavone/flavone core (two benzene rings connected by a 3-carbon chain, one with a carbonyl and an oxygen)
    core_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)C2=Cc3ccccc3O2")
    if not mol.HasSubstructMatch(core_pattern):
        core_pattern = Chem.MolFromSmarts("c1ccccc1C=C2C(=O)c3ccccc3O2") #check flavone core
        if not mol.HasSubstructMatch(core_pattern):
           return False, "No isoflavone/flavone core found"

    # 2. Check for multiple hydroxyl groups (at least 2, often more)
    hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
    if hydroxyl_count < 2:
        return False, f"Too few hydroxyl groups: {hydroxyl_count}"
    
    # 3. Check for glycosylation (at least 1 sugar moiety)
    glycosidic_bond_pattern = Chem.MolFromSmarts("C-O-[C;R]") #check for glycosidic bonds (C-O-C where C is part of ring)
    sugar_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)

    if len(sugar_matches) == 0 and hydroxyl_count<3: # some acrovestones are not glycosides but tend to have 3 or more OH groups
        return False, "No glycosylation detected, and less than 3 hydroxyl groups"

    #additional checks
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 15:
        return False, "Too few carbons for acrovestone"

    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 4:
         return False, "Too few oxygens for acrovestone"

    return True, "Matches acrovestone structural features"