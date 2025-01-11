"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: polymer
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from collections import defaultdict

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is a macromolecule composed of repeating units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - polymers typically have high molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da) to be a polymer"

    # Generate fingerprints to identify repeating units
    # We'll use Morgan fingerprints with radius 2
    fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
    counts = defaultdict(int)
    # Count occurrences of each bit (substructure)
    for bit_id, count in fp.GetNonzeroElements().items():
        counts[bit_id] += count

    # Identify bits that occur multiple times, indicating repeating units
    repeating_units = [bit_id for bit_id, count in counts.items() if count > 2]
    if not repeating_units:
        return False, "No repeating units detected in the molecule"

    # Analyze the degree of polymerization based on repeating units
    total_repeats = sum([counts[bit_id] for bit_id in repeating_units])
    if total_repeats < 10:
        return False, "Not enough repeating units to be considered a polymer"

    return True, "Molecule is a macromolecule with repeating units, indicative of a polymer"