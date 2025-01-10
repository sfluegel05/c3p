"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is an inorganic nutrient essential for the human body.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """

    # Essential mineral elements required by the human body (atomic numbers)
    essential_minerals = {
        11,  # Na - Sodium
        12,  # Mg - Magnesium
        13,  # Al - Aluminum
        15,  # P - Phosphorus
        16,  # S - Sulfur
        17,  # Cl - Chlorine
        19,  # K - Potassium
        20,  # Ca - Calcium
        26,  # Fe - Iron
        29,  # Cu - Copper
        30,  # Zn - Zinc
        53,  # I - Iodine
        55,  # Cs - Caesium
        56,  # Ba - Barium
        57,  # La - Lanthanum
        82,  # Pb - Lead
        24,  # Cr - Chromium
    }

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return False, "Invalid SMILES string"

    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.KekulizeException:
        return False, "Failed to sanitize molecule"

    # Get set of atomic numbers present in the molecule
    atomic_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}

    # Check for presence of essential mineral elements
    minerals_in_mol = atomic_nums & essential_minerals
    if not minerals_in_mol:
        return False, "Does not contain essential mineral elements required by the human body"

    # Check for covalent bonds between essential minerals and carbon
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # If essential mineral is covalently bonded to carbon, exclude
        if (a1.GetAtomicNum() in essential_minerals and a2.GetAtomicNum() == 6) or \
           (a2.GetAtomicNum() in essential_minerals and a1.GetAtomicNum() == 6):
            return False, "Contains covalent bonds between essential minerals and carbon"

    # Count number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Allow small organic counterions (e.g., acetates)
    if num_carbons > 20:
        return False, f"Contains too many carbon atoms ({num_carbons}), likely not a mineral nutrient"

    # Check for aromatic atoms (characteristic of organic molecules)
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > 0:
        return False, "Contains aromatic rings, likely an organic molecule"

    # Passed all checks
    return True, "Molecule is an inorganic compound containing essential mineral nutrients"