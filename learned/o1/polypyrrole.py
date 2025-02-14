"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pyrrole SMARTS pattern
    # Matches a five-membered ring with one nitrogen atom and four carbon atoms
    pyrrole_smarts = '[#6]1:[#6]:[#6]:[#6]:[#7H]:1'  # General pyrrole ring
    pyrrole_pattern = Chem.MolFromSmarts(pyrrole_smarts)

    # Find all matches of the pyrrole substructure
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    pyrrole_count = len(pyrrole_matches)

    # Alternatively, use a substructure query to find five-membered rings with one nitrogen
    pyrrole_count_alt = 0
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    for ring in atom_rings:
        if len(ring) == 5:
            n_count = 0
            c_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    n_count += 1
                    if not atom.GetIsAromatic():
                        continue  # Skip non-aromatic nitrogen
                elif atom.GetAtomicNum() == 6:
                    c_count += 1
                else:
                    break  # Contains other atoms
            if n_count == 1 and c_count == 4:
                pyrrole_count_alt += 1

    # Use the maximum count from both methods
    pyrrole_count = max(pyrrole_count, pyrrole_count_alt)

    if pyrrole_count >= 2:
        return True, f"Contains {pyrrole_count} pyrrole units"
    else:
        return False, f"Contains only {pyrrole_count} pyrrole unit(s), less than 2 required for polypyrrole"