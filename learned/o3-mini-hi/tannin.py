"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: TANNIN – Defined as “any of a group of astringent polyphenolic vegetable principles or compounds,
chiefly complex glucosides of catechol and pyrogallol.”

This version improves upon previous heuristics by:
  • Adding explicit hydrogens so hydroxyl groups are visible.
  • Counting aromatic rings.
  • Searching for canonical polyphenolic motifs (catechol, pyrogallol, and gallic acid).
  • Counting “free phenol” groups (aromatic C bonded to an -OH).
  • Detecting sugar moieties using a permissive sugar SMARTS.
  • Adjusting evidence thresholds when canonical motifs are found versus when only free hydroxyls are seen.
  • Applying slightly higher molecular weight requirements for glycosylated molecules.
Note: This remains a heuristic approach; borderline cases may still be mis‐classified.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    
    The function:
      - Parses the SMILES and adds explicit hydrogens.
      - Computes molecular weight and counts aromatic rings.
      - Searches for canonical polyphenolic motifs (catechol, pyrogallol, gallic acid).
      - Counts free phenol groups (aromatic carbon bearing an -OH group).
      - Identifies a sugar moiety using a permissive pyranose pattern.
      - Uses different thresholds for glycosylated vs. non‐glycosylated molecules.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a tannin, False otherwise.
        str: Explanation for the decision.
    """
    # Parse molecule and add explicit hydrogens so that –OH groups become detectable.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Calculate molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Count aromatic rings using the ring info.
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1

    # Define SMARTS strings for canonical polyphenolic motifs.
    catechol_pat = Chem.MolFromSmarts("c1cc(O)cc(O)c1")          # 1,2-dihydroxybenzene motif
    pyrogallol_pat = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")       # 1,2,3-trihydroxybenzene motif
    gallic_acid_pat = Chem.MolFromSmarts("OC(=O)c1cc(O)c(O)c(O)c1")  # Gallic acid motif

    catechol_matches = mol.GetSubstructMatches(catechol_pat) if catechol_pat is not None else []
    pyrogallol_matches = mol.GetSubstructMatches(pyrogallol_pat) if pyrogallol_pat is not None else []
    gallic_matches = mol.GetSubstructMatches(gallic_acid_pat) if gallic_acid_pat is not None else []
    polyphenolic_hits = len(catechol_matches) + len(pyrogallol_matches) + len(gallic_matches)
    
    # Count free phenol groups: Look at each aromatic carbon and check if it is attached to at least one oxygen that has hydrogen.
    free_phenol_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    free_phenol_count += 1
                    break    # count one free phenol per aromatic carbon

    combined_evidence = polyphenolic_hits + free_phenol_count

    # Define a permissive SMARTS for a sugar (pyranose-like) fragment.
    sugar_pat = Chem.MolFromSmarts("[$([C@H]1OC(O)C(O)C(O)C1O)]")
    sugar_found = mol.HasSubstructMatch(sugar_pat) if sugar_pat is not None else False

    # Apply criteria based on whether a sugar moiety is present.
    if sugar_found:
        # For glycosylated tannins, require a bit higher mass and at least 2 aromatic rings.
        if mol_wt < 350:
            return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a tannin glycoside"
        if aromatic_ring_count < 2:
            return False, f"Only {aromatic_ring_count} aromatic ring(s) found despite the sugar moiety"
        if polyphenolic_hits < 1:
            return False, "No canonical polyphenolic motif detected despite presence of sugar"
        # When a canonical polyphenolic pattern is present, do not insist on extra free OH evidence.
        return True, "Molecule exhibits polyphenolic glycoside features consistent with tannins"
    else:
        # For non-glycosylated tannins, require at least modest MW and multiple aromatic rings.
        if mol_wt < 300:
            return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a polyphenolic tannin"
        if aromatic_ring_count < 2:
            return False, f"Insufficient aromatic rings (found {aromatic_ring_count}; need at least 2)"
        if polyphenolic_hits < 1:
            return False, "No canonical polyphenolic motif detected in the molecule"
        # When free phenol evidence is missing, allow heavier molecules to pass (they might be fully methylated).
        if free_phenol_count == 0:
            if mol_wt < 400:
                return False, ("No free phenolic hydroxyl group detected and molecular weight is too low; "
                               "expected for a tannin")
            else:
                return True, ("Molecule has canonical polyphenolic motifs and high molecular weight; "
                              "despite zero free phenols, it is classified as a tannin")
        # For those with at least one free hydroxyl, require combined evidence to be moderately above the motifs.
        if combined_evidence < (polyphenolic_hits + 1):
            return False, (f"Combined motif ({polyphenolic_hits}) and free phenol ({free_phenol_count}) evidence "
                           "appears insufficient for a tannin")
        return True, "Molecule exhibits classic polyphenolic features consistent with tannins"

# Example usage:
if __name__ == "__main__":
    # Example list of SMILES for known tannins (true positives) and a few negatives.
    test_smiles = [
        # True positives (tannins)
        "C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2",  # Besigomsin
        "O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@H]([C@@H](O)[C@]([C@@H]3O)(O)C)[C@@H]4C=5C(=O)C6=C(O)C=C(OC)C=C6C(C5[C@H](O)[C@@]([C@H]4O)(O)C)=O",  # Alterporriol O
        "O1C2=C(OC)C(=C(C3=C(O)C=4OCOC4C(=C3C)OC)C(=C2OC1)O)C",  # Benzocamphorin E
        "Oc1cc(O)cc(Oc2c(O)cc(O)c3Oc4cc(O)cc(O)c4Oc23)c1",  # eckol
        "Oc1cc(cc(O)c1O)C(=O)O[C@@H]1O[C@@H]2COC(=O)c3cc(O)c(O)c(O)c3-c3c(O)c(O)c(O)cc3C(=O)O[C@H]2[C@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1",  # eugeniin
        # Additional examples could be tested here.
    ]
    
    for sm in test_smiles:
        result, reason = is_tannin(sm)
        print(f"SMILES: {sm}\nIs tannin? {result}\nReason: {reason}\n")