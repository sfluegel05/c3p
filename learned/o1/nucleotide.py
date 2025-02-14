"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:33561 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation
    of the 3' or 5' hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kekulize molecule for accurate aromaticity detection
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except:
        pass

    # Identify sugar rings (ribose or deoxyribose)
    sugar_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@@H](O)[C@H](O[C@H]1CO)O")
    sugars = mol.GetSubstructMatches(sugar_pattern)
    if not sugars:
        # Try deoxyribose pattern
        deoxy_sugar_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@@H](O)[C@H](O[C@H]1CO)")
        sugars = mol.GetSubstructMatches(deoxy_sugar_pattern)
        if not sugars:
            return False, "No sugar ring (ribose or deoxyribose) found"

    if len(sugars) != 1:
        return False, f"Expected one sugar ring, found {len(sugars)}"

    # Identify phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    phosphates = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphates:
        return False, "No phosphate group found"

    # Ensure phosphate is connected to sugar's 3' or 5' oxygen
    phosphate_connected = False
    phosphate_atoms = [atom for match in phosphates for atom in match]
    sugar_ring_atoms = sugars[0]
    sugar_oxygen_atoms = [idx for idx in sugar_ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]

    for o_idx in sugar_oxygen_atoms:
        atom = mol.GetAtomWithIdx(o_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in phosphate_atoms:
                phosphate_connected = True
                break
        if phosphate_connected:
            break

    if not phosphate_connected:
        return False, "Phosphate group not connected to sugar's 3' or 5' oxygen"

    # Identify nitrogenous base connected to sugar
    base_connected = False
    base_atoms = []
    for atom in mol.GetAtomWithIdx(sugar_ring_atoms[0]).GetNeighbors():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() not in sugar_ring_atoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 7 or neighbor.GetAtomicNum() == 6:
                    # Check if neighbor is part of an aromatic ring (base)
                    is_in_ring = neighbor.IsInRing()
                    if is_in_ring:
                        base_connected = True
                        base_atoms.append(neighbor.GetIdx())
                        break
            if base_connected:
                break

    if not base_connected:
        return False, "No nitrogenous base connected to sugar via glycosidic bond"

    # Check that the base is a heterocyclic aromatic ring
    base_query = Chem.MolFromSmarts("[$([n])]:c:c:n:c")  # General pattern for nitrogenous bases
    base_matches = mol.GetSubstructMatches(base_query)

    if not base_matches:
        return False, "Nitrogenous base not found or not aromatic"

    # Ensure only one base is present
    base_rings = [mol.GetAtomWithIdx(idx).GetRingInfo().AtomRings() for idx in base_atoms]
    base_ring_atoms = set()
    for rings in base_rings:
        for ring in rings:
            base_ring_atoms.update(ring)
    if len(base_ring_atoms) > 15:
        return False, "Multiple bases or extended ring systems found"

    # Ensure molecule is not a polymer
    # Count the number of sugars in the molecule
    num_sugars = len(mol.GetSubstructMatches(sugar_pattern)) + len(mol.GetSubstructMatches(deoxy_sugar_pattern))
    if num_sugars > 1:
        return False, "Molecule contains multiple sugar units, possibly an oligonucleotide"

    # Count the number of phosphate groups
    num_phosphates = len(phosphates)
    # Allow for mono-, di-, or triphosphates
    if num_phosphates > 3:
        return False, f"Too many phosphate groups ({num_phosphates}), possibly a polymer"

    return True, "Molecule is a nucleotide"