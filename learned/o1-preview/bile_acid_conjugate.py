"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: Bile Acid Conjugate
"""

from rdkit import Chem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate is a bile acid attached to a hydrophilic group such as glycine, taurine,
    other amino acids, sulfate, glucuronic acid, glucose, other sugars, or coenzyme A.

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

    # Check for steroid nucleus: fused ring system of four rings (6-6-6-5)
    ri = mol.GetRingInfo()
    sssr = ri.AtomRings()
    if len(sssr) < 4:
        return False, "No steroid nucleus found (not enough rings)"

    # Get ring sizes
    ring_sizes = [len(ring) for ring in sssr]

    # Check for ring sizes 6-6-6-5
    six_membered_rings = ring_sizes.count(6)
    five_membered_rings = ring_sizes.count(5)
    if not (six_membered_rings >= 3 and five_membered_rings >= 1):
        return False, "No steroid nucleus found (incorrect ring sizes)"

    # Now check for hydrophilic conjugation groups
    # Define SMARTS patterns for conjugation groups

    # Glycine conjugation (amide bond to glycine)
    glycine_conj_smarts = '[NX3][CX3](=O)[CX4][CX3](=O)[O-,O]'
    glycine_conj = Chem.MolFromSmarts(glycine_conj_smarts)

    # Taurine conjugation (amide bond to taurine)
    taurine_conj_smarts = '[NX3][CX4][CX4][SX4](=O)(=O)[O-,O]'
    taurine_conj = Chem.MolFromSmarts(taurine_conj_smarts)

    # Sulfate conjugation (sulfate ester)
    sulfate_conj_smarts = '[OX2][SX4](=O)(=O)[O-,O]'
    sulfate_conj = Chem.MolFromSmarts(sulfate_conj_smarts)

    # Glucuronic acid conjugation (ester or glycosidic bond to glucuronic acid)
    glucuronic_acid_conj_smarts = 'O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@@H]1C(=O)[O-,O]'
    glucuronic_acid_conj = Chem.MolFromSmarts(glucuronic_acid_conj_smarts)

    # Glucose conjugation (glycosidic bond)
    glucose_conj_smarts = 'O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H]1O'
    glucose_conj = Chem.MolFromSmarts(glucose_conj_smarts)

    # General amino acid conjugation (amide bond to amino acid side chain)
    amino_acid_conj_smarts = '[NX3][CX3](=O)[CX4][CX3](=O)[O-,O]'
    amino_acid_conj = Chem.MolFromSmarts(amino_acid_conj_smarts)

    # List of conjugation patterns
    conjugation_patterns = [
        ('glycine', glycine_conj),
        ('taurine', taurine_conj),
        ('sulfate', sulfate_conj),
        ('glucuronic acid', glucuronic_acid_conj),
        ('glucose', glucose_conj),
        ('amino acid', amino_acid_conj)
    ]

    # Check for conjugation
    for name, pattern in conjugation_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {name} conjugation"

    # Also check for coenzyme A conjugation (simplified pattern)
    coa_conj_smarts = 'NC(=O)CCNC(=O)[CX4][CX3](=O)[O-,O]'
    coa_conj = Chem.MolFromSmarts(coa_conj_smarts)
    if mol.HasSubstructMatch(coa_conj):
        return True, "Contains coenzyme A conjugation"

    return False, "No hydrophilic conjugation found"