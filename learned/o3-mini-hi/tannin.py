"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: TANNIN – Any of a group of astringent polyphenolic vegetable principles or compounds,
chiefly complex glucosides of catechol and pyrogallol.

This updated function attempts to improve the previous heuristic by:
  1. Searching for several polyphenolic motifs (catechol, pyrogallol, and gallic acid).
  2. Counting free phenolic hydroxyl groups (–OH on an aromatic carbon).
  3. Detecting the presence of a sugar moiety.
  4. Using different criteria depending on whether a sugar is present.
  
The idea is that tannins—although diverse—should be rich in aromatic/hydroxyl features.
Some known true tannins (even without sugar) are lower in molecular weight or have fewer aromatic rings,
while many false positives have only “partial” polyphenolic patterns.
Note: This heuristic may still mis‐classify some molecules.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    
    The function uses several criteria:
      1. Molecule parsing and computation of molecular weight.
      2. Counting aromatic rings from ring information.
      3. Searching for polyphenolic substructures:
           - Catechol: 1,2-dihydroxybenzene.
           - Pyrogallol: 1,2,3-trihydroxybenzene.
           - Gallic acid: carboxylated trihydroxybenzene.
      4. Counting free “phenol” groups (an aromatic carbon directly bonded to a hydroxyl group).
      5. Detecting a sugar fragment (a common pyranose-like ring).
      
    Then different thresholds are applied depending on whether a sugar is present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a tannin, False otherwise.
        str: A reason for the classification.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Count aromatic rings using ring info.
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
    
    # Define SMARTS for polyphenolic substructures.
    catechol_pat = Chem.MolFromSmarts("c1cc(O)cc(O)c1")       # 1,2-dihydroxybenzene
    pyrogallol_pat = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")    # 1,2,3-trihydroxybenzene
    gallic_acid_pat = Chem.MolFromSmarts("OC(=O)c1cc(O)c(O)c(O)c1")  # Gallic acid motif

    catechol_matches = mol.GetSubstructMatches(catechol_pat) if catechol_pat is not None else []
    pyrogallol_matches = mol.GetSubstructMatches(pyrogallol_pat) if pyrogallol_pat is not None else []
    gallic_matches = mol.GetSubstructMatches(gallic_acid_pat) if gallic_acid_pat is not None else []
    polyphenolic_hits = len(catechol_matches) + len(pyrogallol_matches) + len(gallic_matches)
    
    # Define SMARTS for free phenol groups.
    # This pattern matches an aromatic carbon directly bonded to an –OH group.
    phenol_pat = Chem.MolFromSmarts("c[OH]")
    phenol_matches = mol.GetSubstructMatches(phenol_pat) if phenol_pat is not None else []
    n_phenols = len(phenol_matches)
    
    # Define a SMARTS for a common sugar (pyranose-like ring with several hydroxyls).
    sugar_pat = Chem.MolFromSmarts("[C@H]1OC(O)C(O)C(O)C1O")
    sugar_found = mol.HasSubstructMatch(sugar_pat) if sugar_pat is not None else False

    # Apply different criteria based on the presence of a sugar moiety.
    if sugar_found:
        # For glycosylated tannins, require a moderate molecular weight,
        # at least one aromatic ring and some polyphenolic/hydroxyl evidence.
        if mol_wt < 300:
            return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a tannin glycoside"
        if aromatic_ring_count < 1:
            return False, f"Insufficient aromatic rings (found {aromatic_ring_count}) even though a sugar was detected"
        if (polyphenolic_hits + n_phenols) < 2:
            return False, "Not enough polyphenolic/hydroxyl evidence despite presence of a sugar moiety"
        return True, "Molecule exhibits polyphenolic glycoside features consistent with tannins"
    else:
        # For non-glycosylated tannins, we try two routes:
        # (a) if one of our classical polyphenolic motifs is found, we allow a lower aromatic ring count,
        # (b) otherwise require at least 2 aromatic rings and sufficient free phenol groups.
        if polyphenolic_hits > 0:
            if mol_wt < 300:
                return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a polyphenolic tannin"
            if aromatic_ring_count < 1:
                return False, f"Insufficient aromatic rings (found {aromatic_ring_count}) for a tannin"
            return True, "Molecule exhibits classic polyphenolic motifs consistent with tannins"
        else:
            # If no direct polyphenolic motif, then require more evidence from overall aromatic rings and free –OH groups.
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
        # True positives:
        "C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2",  # Besigomsin
        "O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@H]([C@@H](O)[C@]([C@@H]3O)(O)C)[C@@H]4C=5C(=O)C6=C(O)C=C(OC)C=C6C(C5[C@H](O)[C@@]([C@H]4O)(O)C)=O",  # Alterporriol O
        "Oc1cc(O)cc(Oc2c(O)cc(O)c3Oc4cc(O)cc(O)c4Oc23)c1",  # eckol (FP previously missed)
        "COCC1=CC(=C(C(=C1C2=C(C(=C(C=C2COC)O)O)O)O)O)O",  # 5-(methoxymethyl)-4-[2,3,4-trihydroxy-6-(methoxymethyl)phenyl]benzene-1,2,3-triol
        "O1C[C@@]23C=4C([C@H](OC(=O)/C(/C)=C\\C)[C@@H]([C@@H](CC2=CC(OC)=C(OC)C3=O)C)C)=CC=5OCOC5C14",  # Heteroclitin D
        
        # False positives (expected to be rejected):
        "O=C(OC1=C(C(O)=C(C(=O)O)C(=C1C)C)C)C2=C(OC)C(=C(OC(=O)C3=C(O)C=C(O)C=C3C)C=C2C)C",  # Thielavin Z5
        "Oc1ccc(cc1OS([O-])(=O)=O)-c1oc2cc([O-])cc(O)c2c(=O)c1OS([O-])(=O)=O"  # quercetin 3,3'-bissulfate(3-)
    ]
    
    for sm in test_smiles:
        result, reason = is_tannin(sm)
        print(f"SMILES: {sm}\nIs tannin? {result}\nReason: {reason}\n")