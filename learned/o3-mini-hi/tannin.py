"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: TANNIN – Any of a group of astringent polyphenolic vegetable principles or compounds,
chiefly complex glucosides of catechol and pyrogallol.

This improved function uses the following approach:
  1. Parses the input SMILES.
  2. Computes molecular weight and counts aromatic rings.
  3. Searches for specific polyphenolic motifs (catechol, pyrogallol, gallic acid).
  4. Counts “free” phenolic hydroxyl groups (an aromatic carbon directly attached to –OH).
  5. Detects a sugar moiety by matching a common pyranose-like fragment.
  
For glycosylated compounds the overall evidence threshold is set to require slightly higher
combined polyphenol evidence (motif hits + free phenols). For non-glycosylated compounds, if none
of the canonical polyphenolic motifs are directly found, we require at least two aromatic rings and
three free phenolic –OH groups.

While still heuristic, these criteria aim to reduce false positives (molecules with only partial
polyphenolic motifs) while correctly classifying known tannins.

Note: This approach may still miss some borderline cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    
    The following checks are applied:
      - Valid SMILES and computed molecular weight.
      - Count of aromatic rings.
      - Identification of common polyphenolic motifs:
           • Catechol: 1,2-dihydroxybenzene.
           • Pyrogallol: 1,2,3-trihydroxybenzene.
           • Gallic acid: carboxylated trihydroxybenzene.
      - Count of free phenol groups (aromatic carbon with an –OH).
      - Detection of a sugar moiety (via a pyranose-like fragment).
    
    Different thresholds are applied when a sugar is present versus when it is not.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a tannin, False otherwise.
        str: A reason for the classification.
    """
    # Parse molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Compute molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    # Count aromatic rings using the molecule's ring info.
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1

    # Define SMARTS for polyphenolic motifs.
    catechol_pat = Chem.MolFromSmarts("c1cc(O)cc(O)c1")         # 1,2-dihydroxybenzene
    pyrogallol_pat = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")      # 1,2,3-trihydroxybenzene
    gallic_acid_pat = Chem.MolFromSmarts("OC(=O)c1cc(O)c(O)c(O)c1") # Gallic acid motif

    catechol_matches = mol.GetSubstructMatches(catechol_pat) if catechol_pat is not None else []
    pyrogallol_matches = mol.GetSubstructMatches(pyrogallol_pat) if pyrogallol_pat is not None else []
    gallic_matches = mol.GetSubstructMatches(gallic_acid_pat) if gallic_acid_pat is not None else []
    polyphenolic_hits = len(catechol_matches) + len(pyrogallol_matches) + len(gallic_matches)
    
    # Define a SMARTS for free phenol groups: aromatic carbon bonded to hydroxyl.
    phenol_pat = Chem.MolFromSmarts("c[OH]")
    phenol_matches = mol.GetSubstructMatches(phenol_pat) if phenol_pat is not None else []
    n_phenols = len(phenol_matches)
    
    # Define a SMARTS for a common sugar fragment (pyranose ring with several hydroxyls).
    sugar_pat = Chem.MolFromSmarts("[C@H]1OC(O)C(O)C(O)C1O")
    sugar_found = mol.HasSubstructMatch(sugar_pat) if sugar_pat is not None else False

    # Apply different criteria based on the presence of a sugar.
    if sugar_found:
        # Glycosylated tannins:
        if mol_wt < 300:
            return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a tannin glycoside"
        if aromatic_ring_count < 1:
            return False, f"Insufficient aromatic rings (found {aromatic_ring_count}) despite sugar moiety"
        # Increase the evidence threshold: require at least 3 combined hits.
        if (polyphenolic_hits + n_phenols) < 3:
            return False, "Not enough polyphenolic/hydroxyl evidence despite presence of a sugar moiety"
        return True, "Molecule exhibits polyphenolic glycoside features consistent with tannins"
    else:
        # Non-glycosylated tannins:
        if polyphenolic_hits > 0:
            # At least one classical polyphenolic motif was found.
            if mol_wt < 300:
                return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a polyphenolic tannin"
            if aromatic_ring_count < 1:
                return False, f"Insufficient aromatic rings (found {aromatic_ring_count}) for a tannin"
            return True, "Molecule exhibits classic polyphenolic motifs consistent with tannins"
        else:
            # If none of the clear motifs are present, require more evidence.
            if mol_wt < 300:
                return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a tannin lacking sugar"
            if aromatic_ring_count < 2:
                return False, f"Insufficient aromatic rings (found {aromatic_ring_count}, need at least 2) for a tannin without sugar"
            if n_phenols < 3:
                return False, f"Too few free phenolic hydroxyl groups (found {n_phenols}); expected at least 3 for a tannin without sugar"
            return True, "Molecule exhibits rich polyphenolic features consistent with tannins"

# Example usage:
if __name__ == "__main__":
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
        "S(OC=1C(O)=CC(C(OC2C(OC(=O)C3=CC(O)=C(O)C(O)=C3)C(OC(O)C2O)CO)=O)=CC1O)(O)(=O)=O",  # [4-({[2,3-dihydroxy-6-(hydroxymethyl)-5-(3,4,5-trihydroxybenzoyloxy)oxan-4-yl]oxy}carbonyl)-2,6-dihydroxyphenyl]oxidanesulfonic acid
        "O[C@@H]1([C@@H](O)[C@H](OC(=O)c1cc(O)c(O)c(O)c1)C)[C@H](O)OC1=C(O)c(O)c(O)c(O)C=C1",  # Thonningianin A (simplified example)
        "O[C@H]1Cc2c(O)cc(O)c([C@H]3[C@H](O)[C@H](Oc4cc(O)c(O)c34)c3cc(O)c(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1",  # (+)-gallocatechin-(4alpha->8)-(+)-catechin
        "Oc1cc(cc(O)c1O)C(=O)O[C@@H]1O[C@@H]2COC(=O)c3cc(O)c(O)c(O)c3-c3c(O)c(O)c(O)cc3C(=O)O[C@H]2[C@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1",  # Isoterchebin
        "COc1cc2[C@H](OC(=O)c3ccccc3)[C@@](C)(O)[C@@H](C)Cc3cc4OCOc4c(OC)c3-c2c(OC)c1OC",  # Schisantherin A
        "O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@H](O)[C@]([C@H]3O)(O)C)[C@H]4C=5C(=O)C6=C(O)C=C(OC)C=C6C(C5[C@H](O)[C@@]([C@H]4O)(O)C)=O",  # Alterporriol S
        "O[C@@H]1Cc2c(O)cc(O)c([C@@H]3[C@@H](O)[C@H](Oc4cc(O)c(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1",  # procyanidin B2
        "O1C[C@@]23C=4C([C@H](OC(=O)/C(/C)=C\\C)[C@@H]([C@@H](CC2=CC(OC)=C(OC)C3=O)C)C)=CC=5OCOC5C14",  # Heteroclitin D
        
        # False positives (expected to be rejected):
        "O=C(OC1=C(C(O)=C(C(=O)O)C(=C1C)C)C)C2=C(OC)C(=C(OC(=O)C3=C(O)C=C(O)C=C3C)C=C2C)C",  # Thielavin Z5
        "Oc1ccc(cc1OS([O-])(=O)=O)-c1oc2cc([O-])cc(O)c2c(=O)c1OS([O-])(=O)=O"  # quercetin 3,3'-bissulfate(3-)
    ]
    
    for sm in test_smiles:
        result, reason = is_tannin(sm)
        print(f"SMILES: {sm}\nIs tannin? {result}\nReason: {reason}\n")