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

    # 1. Check for a flavone or isoflavone-like core
    # Define core substructures
    # Benzene and pyrone
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    pyrone_pattern = Chem.MolFromSmarts("c1cc(=O)oc1") #or "c1cc(=O)cc1"
    if not (mol.HasSubstructMatch(benzene_pattern) and mol.HasSubstructMatch(pyrone_pattern)):
          return False, "No flavone/isoflavone-like core found (no benzene-pyrone combination)"

    # Check if the benzene is connected directly to the pyrone, A and B ring structure
    a_ring_pattern = Chem.MolFromSmarts("c1cc(cc(c1)O)C(=O)c2ccccc2")  # A ring with OH group.
    b_ring_pattern = Chem.MolFromSmarts("c1ccccc1O") #B-ring with OH group, this is very simple for now.
    if not (mol.HasSubstructMatch(a_ring_pattern)): #require at least an A ring, and one OH on this ring.
         return False, "A ring is missing"

    # 2. Check for glycosylation (optional)
    # Look for glycosidic bonds, C-O-C between a sugar and a core carbon
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C;R]-O-[C;R][OX2]")
    sugar_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    is_glycosylated = len(sugar_matches) > 0

    # 3. Check for multiple hydroxyl groups on rings
    # count hydroxyls on ring
    hydroxyl_pattern = Chem.MolFromSmarts("c1ccccc1O")
    ring_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if ring_hydroxyls < 1: # At least one hydroxyl is required on a ring
        return False, "Must have at least 1 hydroxyl on ring."

    # 4. Molecular Weight check, acrovestones tend to be heavier
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for acrovestone"

    #Final check
    if ring_hydroxyls > 0 or is_glycosylated:
         return True, "Matches acrovestone structural features"


    return False, "Does not match acrovestone characteristics"