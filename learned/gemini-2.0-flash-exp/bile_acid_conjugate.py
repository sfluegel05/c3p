"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    Bile acid conjugates have a steroidal core and are conjugated to a hydrophilic/charged group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More generic bile acid core pattern
    # Allows for variations at ring junctions, double bonds and substitutions
    bile_acid_core_pattern = Chem.MolFromSmarts("[C]1[C][C][C]2[C]3[C]([C]([C]1[C]2)([C])([C])[C])[C][C][C]4[C]([C]3[C])[C][C][C]4")
    
    if not mol.HasSubstructMatch(bile_acid_core_pattern):
        return False, "No bile acid core structure detected."

    # Define common conjugation groups using SMARTS patterns

    # Broad amino acid pattern (including glycine, taurine etc)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])")
    
    # Taurine specifically (amino sulfonic acid):
    taurine_pattern = Chem.MolFromSmarts("NCCS(=O)(=O)[O,H]")
    
    # Sulfates
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O,H]")

    # Glucuronides
    glucuronide_pattern = Chem.MolFromSmarts("C1[C@H]([C@H]([C@@H]([C@H](O1)O)O)O)C(=O)O")
   
   # Check for the presence of at least one conjugation group
    if not (mol.HasSubstructMatch(amino_acid_pattern) or \
            mol.HasSubstructMatch(taurine_pattern) or \
            mol.HasSubstructMatch(sulfate_pattern) or \
            mol.HasSubstructMatch(glucuronide_pattern)):
        return False, "No common conjugation group detected."
    
    # Molecular Weight filter - bile acid conjugates are generally large.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for bile acid conjugate"

    # Sanity check for number of carbon and oxygen
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for bile acid conjugate."

    if o_count < 3:
        return False, "Too few oxygens for bile acid conjugate."


    return True, "Contains a bile acid core and at least one common conjugation group."