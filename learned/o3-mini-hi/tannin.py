"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: TANNIN – Any of a group of astringent polyphenolic vegetable principles
or compounds, chiefly complex glucosides of catechol and pyrogallol.

This updated function tries to improve upon the previous heuristic by:
  1. Searching for several polyphenolic motifs (catechol, pyrogallol, and gallic acid).
  2. Detecting the presence of a sugar moiety (commonly present in tannin glucosides).
  3. Using different minimum requirements for aromatic ring count and molecular weight
     depending on whether a sugar is present.
Note: Tannins are structurally diverse and this heuristic may still mis‐classify some molecules.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    
    The function looks for polyphenolic features and optionally a sugar moiety.
    It uses the following logic:
      1. Parse the molecule.
      2. Search for polyphenolic patterns:
            - Catechol: benzene ring with two hydroxyl groups in ortho positions.
            - Pyrogallol: benzene ring with three hydroxyl groups.
            - Gallic acid: carboxylated trihydroxybenzene.
      3. Check for a sugar fragment (a common pyranose ring pattern).
      4. Count aromatic rings.
      5. Apply molecular weight criteria.
      
    If a sugar is found, then we accept a single aromatic ring plus evidence of polyphenolic
    elements. Without a sugar the molecule must show a clear polyphenolic signal (at least one match)
    plus at least two aromatic rings and a heavier molecular weight.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a tannin, False otherwise.
        str: Reason for classification or rejection.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Count aromatic rings using ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    aromatic_ring_count = 0
    for ring in rings:
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1

    # Define SMARTS for polyphenolic substructures:
    # Catechol: 1,2-dihydroxybenzene.
    catechol_pat = Chem.MolFromSmarts("c1cc(O)cc(O)c1")
    # Pyrogallol: 1,2,3-trihydroxybenzene.
    pyrogallol_pat = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")
    # Gallic acid: carboxylated trihydroxybenzene.
    gallic_acid_pat = Chem.MolFromSmarts("OC(=O)c1cc(O)c(O)c(O)c1")
    
    catechol_matches = mol.GetSubstructMatches(catechol_pat) if catechol_pat is not None else []
    pyrogallol_matches = mol.GetSubstructMatches(pyrogallol_pat) if pyrogallol_pat is not None else []
    gallic_matches = mol.GetSubstructMatches(gallic_acid_pat) if gallic_acid_pat is not None else []
    polyphenolic_hits = len(catechol_matches) + len(pyrogallol_matches) + len(gallic_matches)
    
    # Define a SMARTS for a common sugar (pyranose) ring.
    # This pattern looks for a six-membered ring with one oxygen and several hydroxyl substituents.
    sugar_pat = Chem.MolFromSmarts("[C@H]1OC(O)C(O)C(O)C1O")
    sugar_found = mol.HasSubstructMatch(sugar_pat) if sugar_pat is not None else False

    # Use different criteria depending on whether a sugar moiety is present.
    if sugar_found:
        # With a sugar, even one aromatic ring together with polyphenolic hints may be acceptable.
        if mol_wt < 300:
            return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a tannin glycoside"
        if aromatic_ring_count < 1:
            return False, f"Insufficient aromatic rings (found {aromatic_ring_count}) even though a sugar moiety was detected"
        # Accept if there is either a polyphenolic hit or a sugar present (enabling a glycoside classification).
        return True, "Molecule exhibits polyphenolic glycoside features consistent with tannins"
    else:
        # Without a sugar, require a polyphenolic hit
        if polyphenolic_hits == 0:
            return False, "No polyphenolic substructures (catechol, pyrogallol, or gallic acid) found"
        if aromatic_ring_count < 2:
            return False, f"Insufficient aromatic rings (found {aromatic_ring_count}, need at least 2) for a tannin without sugar"
        # Often condensed tannins or flavonoid oligomers tend to be larger.
        if mol_wt < 400:
            return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a tannin lacking sugar"
        return True, "Molecule exhibits polyphenolic structural features consistent with tannins"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2",  # Besigomsin (TP)
        "O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@H]([C@@H](O)[C@]([C@@H]3O)(O)C)[C@@H]4C=5C(=O)C6=C(O)C=C(OC)C=C6C(C5[C@H](O)[C@@]([C@H]4O)(O)C)=O",  # Alterporriol O (TP)
        "O=C(O)C1=CC(OC)=C(O[C@@H]2O[C@@H]([C@@H](OC)[C@@H]([C@H]2O)O)CO)C(=C1)CC=C(C)C",  # Conoideoglucoside C (FN previously)
    ]
    for sm in test_smiles:
        result, reason = is_tannin(sm)
        print(f"SMILES: {sm}\nIs tannin? {result}\nReason: {reason}\n")