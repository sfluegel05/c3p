"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: TANNIN – Any of a group of astringent polyphenolic vegetable principles or compounds,
chiefly complex glucosides of catechol and pyrogallol.

This function attempts to improve on previous heuristics by:
  • Explicitly adding H atoms so that hydroxyl groups on aromatic rings are identified.
  • Counting aromatic rings.
  • Searching for canonical polyphenolic motifs (catechol, pyrogallol, gallic acid).
  • Counting “free phenol” groups by checking aromatic carbons bonded to hydroxyl groups.
  • Detecting the presence of a sugar moiety using a more permissive pyranose-like SMARTS.
  • Applying different evidence thresholds for glycosylated vs. non‐glycosylated molecules.
Note: This remains a heuristic approach and borderline cases may still be mis‐classified.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.

    The function:
      - Parses the SMILES and adds explicit hydrogens.
      - Computes molecular weight and counts aromatic rings.
      - Searches for polyphenolic motifs (catechol, pyrogallol, gallic acid).
      - Counts free phenol groups (aromatic carbon directly bonded to a hydroxyl).
      - Identifies a sugar moiety using a pyranose-like fragment SMARTS.
      - Applies different thresholds for glycosylated versus non-glycosylated molecules.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if classified as a tannin, False otherwise.
        str: Explanation for the decision.
    """
    # Parse molecule and add H atoms (so we see –OH groups)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Compute molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Count aromatic rings using the ring info
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1

    # Define SMARTS for canonical polyphenolic motifs.
    catechol_pat = Chem.MolFromSmarts("c1cc(O)cc(O)c1")         # 1,2-dihydroxybenzene
    pyrogallol_pat = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")      # 1,2,3-trihydroxybenzene
    gallic_acid_pat = Chem.MolFromSmarts("OC(=O)c1cc(O)c(O)c(O)c1") # Gallic acid motif
    
    catechol_matches = mol.GetSubstructMatches(catechol_pat) if catechol_pat is not None else []
    pyrogallol_matches = mol.GetSubstructMatches(pyrogallol_pat) if pyrogallol_pat is not None else []
    gallic_matches = mol.GetSubstructMatches(gallic_acid_pat) if gallic_acid_pat is not None else []
    polyphenolic_hits = len(catechol_matches) + len(pyrogallol_matches) + len(gallic_matches)
    
    # Count free phenol groups – for each aromatic carbon, if it has at least one oxygen neighbor that is part of an -OH, count it.
    free_phenol_count = 0
    for atom in mol.GetAtoms():
        # Check if the atom is an aromatic carbon.
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            # Look for an oxygen neighbor that is part of an -OH (has at least one hydrogen)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen has at least one hydrogen attached.
                    if nbr.GetTotalNumHs() > 0:
                        free_phenol_count += 1
                        break  # count at most one free phenol per aromatic carbon

    # Define a SMARTS for a common sugar fragment (pyranose ring with several hydroxyls).
    # This SMARTS is made a bit more permissive.
    sugar_pat = Chem.MolFromSmarts("[$([C@H]1OC(O)C(O)C(O)C1O)]")
    sugar_found = mol.HasSubstructMatch(sugar_pat) if sugar_pat is not None else False

    # Compute a combined polyphenol evidence score.
    combined_evidence = polyphenolic_hits + free_phenol_count

    # Now apply thresholds.
    if sugar_found:
        # For glycosylated tannins:
        if mol_wt < 300:
            return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a tannin glycoside"
        if aromatic_ring_count < 1:
            return False, f"Insufficient aromatic rings (found {aromatic_ring_count}) despite sugar moiety"
        # Lower threshold due to sugar: require at least 2 combined motif/free phenol hits.
        if combined_evidence < 2:
            return False, "Not enough polyphenolic/hydroxyl evidence despite presence of a sugar moiety"
        return True, "Molecule exhibits polyphenolic glycoside features consistent with tannins"
    else:
        # For non-glycosylated tannins:
        if mol_wt < 300:
            return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a polyphenolic tannin"
        if aromatic_ring_count < 2:
            return False, f"Insufficient aromatic rings (found {aromatic_ring_count}, need at least 2) for a tannin without sugar"
        # If a canonical motif is present, require that free phenol evidence together with motifs is moderately high.
        if polyphenolic_hits > 0:
            if combined_evidence < (polyphenolic_hits + 2):
                return False, (f"Combined motif ({polyphenolic_hits}) and free phenol ({free_phenol_count}) evidence "
                               "are insufficient despite detection of a polyphenolic motif")
        else:
            # If no canonical motif then require at least 3 free evidence hits.
            if combined_evidence < 3:
                return False, f"Too few free phenolic hydroxyl groups (found {free_phenol_count}); expected at least 3 for a tannin lacking sugar"
        return True, "Molecule exhibits classic polyphenolic motifs consistent with tannins"

# Example usage:
if __name__ == "__main__":
    # Some example SMILES from known tannins and some expected negatives:
    test_smiles = [
        # True positives (tannins)
        "C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2",  # Besigomsin
        "O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@H]([C@@H](O)[C@]([C@@H]3O)(O)C)[C@@H]4C=5C(=O)C6=C(O)C=C(OC)C=C6C(C5[C@H](O)[C@@]([C@H]4O)(O)C)=O",  # Alterporriol O
        "O1C2=C(OC)C(=C(C3=C(O)C=4OCOC4C(=C3C)OC)C(=C2OC1)O)C",  # Benzocamphorin E
        "Oc1cc(O)cc(Oc2c(O)cc(O)c3Oc4cc(O)cc(O)c4Oc23)c1",  # eckol
        "Oc1cc(cc(O)c1O)C(=O)O[C@@H]1O[C@@H]2COC(=O)c3cc(O)c(O)c(O)c3-c3c(O)c(O)c(O)cc3C(=O)O[C@H]2[C@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1",  # eugeniin
        "Oc1cc(O)cc(Oc2c(O)cc(O)c3Oc4cc(Oc5c(O)cc(Oc6c(O)cc(O)c7Oc8cc(O)cc(O)c8Oc67)cc5O)cc(O)c4Oc23)c1",  # dieckol
        "COCC1=CC(=C(C(=C1C2=C(C(=C(C=C2COC)O)O)O)O)O)O",  # 5-(methoxymethyl)-4-[2,3,4-trihydroxy-6-(methoxymethyl)phenyl]benzene-1,2,3-triol
        "COc1c(O)cc(cc1O)[C@H]1Oc2c(C[C@@H]1O)c(O)cc(O)c2[C@@H]1[C@@H](O)[C@H](Oc2cc(O)cc(O)c12)c1ccc(O)c(O)c1",  # epicatechin-(4beta->8)-4'-O-methylgallocatechin
        "O[C@@H]1[C@H](Oc2cc(O)cc(O)c2[C@H]1c1c(O)cc(O)c2C[C@@H](OC(=O)c3cc(O)c(O)c(O)c3)[C@H](Oc12)c1ccc(O)c(O)c1)c1ccc(O)c(O)c1",  # procyanidin B4 3'-O-gallate
        "O[C@H]1Cc2c(O)cc(O)c([C@H]3[C@H](O)[C@H](Oc4c([C@H]5[C@H](O)[C@H](Oc6cc(O)cc(O)c56)c5ccc(O)c(O)c5)c(O)cc(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1",  # procyanidin C2
        "O=C(OC1=CC(O)=CC(=C1)CCCCC)C2=C(O[C@@H]3O[C@H](C(=O)O)[C@H](O)[C@@H]([C@H]3O)O)C=C(O)C=C2CCCCC",  # Ascotricin B
        "COc1cc2CC(C)C(C)Cc3cc(OC)c(OC)c(OC)c3-c2c(OC)c1OC",  # Deoxyschizandrin
        "S(OC=1C(O)=CC(C(OC2C(OC(=O)C3=CC(O)=C(O)C(O)=C3)C(OC(O)C2O)CO)=O)=CC1O)(O)(=O)=O",  # [4-({...)]oxidanesulfonic acid
        # False positives (expected to be rejected)
        "O=C(OC1=C(C(O)=C(C(=O)O)C(=C1C)C)C)C2=C(OC)C(=C(OC(=O)C3=C(O)C=C(O)C=C3C)C=C2C)C",  # Thielavin Z5
        "Oc1ccc(cc1OS([O-])(=O)=O)-c1oc2cc([O-])cc(O)c2c(=O)c1OS([O-])(=O)=O"  # quercetin 3,3'-bissulfate(3-)
    ]
    
    for sm in test_smiles:
        result, reason = is_tannin(sm)
        print(f"SMILES: {sm}\nIs tannin? {result}\nReason: {reason}\n")