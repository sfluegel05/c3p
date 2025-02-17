"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: TANNIN – Any of a group of astringent polyphenolic vegetable principles
or compounds, chiefly complex glucosides of catechol and pyrogallol.

This function attempts to identify tannins by:
  1. Checking for polyphenolic substructures (catechol or pyrogallol motifs).
  2. Verifying that the molecule has at least two aromatic rings.
  3. Checking that the molecular weight is in a typical range for tannins (here > 300 Da).

Note: This is a heuristic approach. The chemistry of tannins is complex;
thus, the program may not be perfect for every structure.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    
    The function looks for polyphenolic features:
      - Presence of catechol (1,2-dihydroxybenzene) or pyrogallol (1,2,3-trihydroxybenzene) motifs.
      - At least two aromatic rings.
      - A molecular weight above a minimal cutoff (here chosen as 300 Da) to avoid very small polyphenols.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a tannin, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define substructure patterns for catechol and pyrogallol.
    # Catechol: benzene ring with two hydroxyl groups in ortho positions.
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)cc(O)c1")
    # Pyrogallol: benzene ring with three hydroxyl groups.
    pyrogallol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")
    
    matches_catechol = mol.GetSubstructMatches(catechol_pattern) if catechol_pattern is not None else []
    matches_pyrogallol = mol.GetSubstructMatches(pyrogallol_pattern) if pyrogallol_pattern is not None else []
    
    if not matches_catechol and not matches_pyrogallol:
        return False, "No catechol/pyrogallol substructures found"
    
    # Count aromatic rings in the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    aromatic_ring_count = 0
    for ring in rings:
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
    if aromatic_ring_count < 2:
        return False, f"Insufficient aromatic rings (found {aromatic_ring_count}, need at least 2)"
    
    # Check a minimal molecular weight cutoff.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a typical tannin"
    
    return True, "Molecule exhibits polyphenolic structural features consistent with tannins"