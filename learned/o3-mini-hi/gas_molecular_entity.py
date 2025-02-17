"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: Gas Molecular Entity
Definition: Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa)
Note: This classifier uses heuristics based on element composition, molecular weight, and the overall size of the molecule.
It may not be 100% accurate.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a given molecule (specified as a SMILES string) likely corresponds to a gas molecular entity at STP,
    based on its composition and size.

    Heuristic criteria used:
      1. The molecule must be parseable.
      2. It must consist only of main group elements (allowed atomic numbers: H, He, B, C, N, O, F, Ne, Na, Mg,
         Al, Si, P, S, Cl, Ar, K, Ca, Br, Kr, I, Xe, Rn).
      3. Its molecular weight must be below a threshold (here, 300 Da) since most gaseous molecules are small.
      4. The molecule should not be too large (here, we require no more than 25 heavy atoms).
      5. The molecule must be contiguous (a single molecular entity).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a gas molecular entity, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Allowed main group elements (using their atomic numbers):
    # H (1), He (2), B (5), C (6), N (7), O (8), F (9), Ne (10),
    # Na (11), Mg (12), Al (13), Si (14), P (15), S (16), Cl (17),
    # Ar (18), K (19), Ca (20), Br (35), Kr (36), I (53), Xe (54), Rn (86)
    allowed_atomic_numbers = {1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                              14, 15, 16, 17, 18, 19, 20, 35, 36, 53, 54, 86}

    # Check that every atom is a main group element
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_numbers:
            return False, f"Atom {atom.GetSymbol()} (atomic number {atom.GetAtomicNum()}) is not a main group element"

    # Calculate the molecular weight using an exact method
    mw = Descriptors.ExactMolWt(mol)
    # Most small gas molecules have low molecular weight; we use 300 Da as a (heuristic) upper threshold
    if mw >= 300:
        return False, f"Molecular weight {mw:.1f} Da is too high for a typical gas molecular entity"

    # Count heavy atoms (atoms other than hydrogen)
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if len(heavy_atoms) > 25:
        return False, f"Too many heavy atoms ({len(heavy_atoms)}) for a typical gas molecular entity"

    # Ensure the molecule is a single connected species (not a salt or a mixture)
    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        return False, "Molecule is disconnected; it should be a single molecular entity to be considered a gas"

    # The molecule passed all the heuristic criteria
    return True, f"Molecule (MW {mw:.1f} Da, {len(heavy_atoms)} heavy atoms) is classified as a gas molecular entity"