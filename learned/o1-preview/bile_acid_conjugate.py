"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: Bile Acid Conjugate
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate is a bile acid attached to a functional group that adds hydrophilicity or charge,
    such as glycine, taurine, other amino acids, sulfate, glucuronic acid, glucose, or coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for bile acid core structure
    # Bile acids are steroids with a specific ring system and side chain
    bile_acid_core_smarts = "C1[C@H]2CC[C@H]3C[C@@H](O)[C@]4(C)CC[C@@H](O)CC4=C[C@H]3[C@@]2(C)CC[C@H]1C(=O)O"
    bile_acid_core = Chem.MolFromSmarts(bile_acid_core_smarts)
    if not mol.HasSubstructMatch(bile_acid_core):
        return False, "No bile acid core structure found"

    # Check for conjugation groups attached to the side chain
    # Glycine conjugation (amide bond to glycine)
    glycine_conj_smarts = "[NX3][C@@H](C(=O)[O-,O])[C@@H](C(=O)[O-,O])N"
    glycine_conj = Chem.MolFromSmarts(glycine_conj_smarts)

    # Taurine conjugation (amide bond to taurine)
    taurine_conj_smarts = "[NX3][C@@H](C(=O)[O-,O])[C@@H](CCS(=O)(=O)[O-,O])[N]"
    taurine_conj = Chem.MolFromSmarts(taurine_conj_smarts)

    # Sulfate conjugation attached to hydroxyl groups
    sulfate_conj_smarts = "[OX2H][SX4](=O)(=O)[O-]"
    sulfate_conj = Chem.MolFromSmarts(sulfate_conj_smarts)

    # Glucuronic acid conjugation (glycosidic bond to glucuronic acid)
    glucuronic_acid_conj_smarts = "O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@@H]1C(=O)[O-,O]"
    glucuronic_acid_conj = Chem.MolFromSmarts(glucuronic_acid_conj_smarts)

    # Glucose conjugation (glycosidic bond)
    glucose_conj_smarts = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"
    glucose_conj = Chem.MolFromSmarts(glucose_conj_smarts)

    # General amino acid conjugation (amide bond to amino acid side chain)
    amino_acid_conj_smarts = "[NX3][C@@H](C(=O)[O-,O])[C@H](O)[C@H](O)[C@@H](C(=O)[O-,O])[N]"
    amino_acid_conj = Chem.MolFromSmarts(amino_acid_conj_smarts)

    # Coenzyme A conjugation (thioester bond to coenzyme A)
    coa_conj_smarts = "C(=O)SCCNC(=O)[CX4]"
    coa_conj = Chem.MolFromSmarts(coa_conj_smarts)

    # List of conjugation patterns
    conjugation_patterns = [
        ('glycine', glycine_conj),
        ('taurine', taurine_conj),
        ('sulfate', sulfate_conj),
        ('glucuronic acid', glucuronic_acid_conj),
        ('glucose', glucose_conj),
        ('amino acid', amino_acid_conj),
        ('coenzyme A', coa_conj)
    ]

    # Check for conjugation
    for name, pattern in conjugation_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {name} conjugation"

    # Also check for other amino acids conjugated via amide bond
    # Pattern for amide bond at side chain
    amide_conj_smarts = "[CX3](=O)[NX3][CX4][CX3](=O)[O-,O]"
    amide_conj = Chem.MolFromSmarts(amide_conj_smarts)
    if mol.HasSubstructMatch(amide_conj):
        return True, "Contains amino acid conjugation"

    return False, "No bile acid conjugation found"