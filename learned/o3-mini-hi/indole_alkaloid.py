"""
Classifies: CHEBI:38958 indole alkaloid
"""
#!/usr/bin/env python
"""
Classifies: Indole Alkaloid – an alkaloid that contains an indole skeleton.
An indole core is defined as the bicyclic system consisting of a benzene ring fused to a pyrrole.
In our approach we (a) remove stereochemical annotations to improve substructure matching,
(b) use two indole SMARTS patterns to catch both free‐NH and N‐substituted indoles,
and (c) apply a peptide filter (by counting simple amide bonds vs. number of rings).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.

    The function checks that:
      1. The SMILES string is valid and the molecule can be sanitized.
      2. The stereochemical information is removed (so as not to miss matches).
      3. The molecule contains a fused bicyclic indole core. Two different SMARTS patterns 
         (one for a classic indole with free NH and one for an N‐substituted indole) are used.
      4. A simple peptide filter is applied: if the molecule contains two or more amide bonds 
         (i.e. “C(=O)N”) and has fewer than 4 rings overall, it is more likely a peptide side‐chain 
         (e.g. from tryptophan in a peptide) than a genuine indole alkaloid.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule qualifies as an indole alkaloid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Try sanitizing the molecule.
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return False, "Molecule could not be sanitized"
    
    # Remove stereochemistry: generate a new molecule with isomeric information stripped.
    # This improves matching when the only differences are chiral annotations.
    mol_stripped = Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=False))
    if mol_stripped is None:
        return False, "Failed to remove stereochemistry for matching"
    
    # --- Peptide filter ---
    # Many peptides (or small peptide fragments) contain one indole side‐chain (from tryptophan)
    # but are not classified as indole alkaloids.
    # Here we check for simple amide bonds via the SMARTS "C(=O)N".
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amide = len(amide_matches)
    # Count the number of rings in the original (with stereochem) molecule.
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    # Heuristically, if there are two or more amide bonds and fewer than four rings overall,
    # we assume the molecule is peptide-like.
    if num_amide >= 2 and num_rings < 4:
        return False, ("Molecule appears to be peptide-like (has {} amide bonds and only {} rings) "
                       "and is not considered an indole alkaloid.".format(num_amide, num_rings))
    
    # --- Indole SMARTS patterns ---
    # Pattern 1: classic indole ring (benzene fused to a pyrrole with free NH)
    indole_smarts1 = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
    # Pattern 2: indole ring with nitrogen substituted (no explicit H on N)
    indole_smarts2 = Chem.MolFromSmarts("c1ccc2c(c1)[n]ccc2")
    
    # Check if either pattern is matched in the stereochem-stripped molecule.
    if mol_stripped.HasSubstructMatch(indole_smarts1) or mol_stripped.HasSubstructMatch(indole_smarts2):
        return True, ("Molecule contains an indole core (benzene fused to a pyrrole ring) "
                      "and can be classified as an indole alkaloid.")
    
    return False, "Molecule does not contain an indole core characteristic of indole alkaloids."

# Example usage:
if __name__ == "__main__":
    # Examples (both true positives and known false positives/negatives)
    test_examples = {
        "staurosporine": "CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)N1C3=C(C=CC=C3)C3=C1C1=C(C4=C(C=CC=C4)N21)C1=C3CNC1=O",
        "3alpha(S)-strictosidine": "[H][C@@]1(C[C@H]2[C@@H](C=C)[C@H](O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)OC=C2C(=O)OC)NCCC2=C1NC1=C2C=CC=C1",
        "folicanthine": "CN1CCC2(C1N(C)c1ccccc21)C12CCN(C)C1N(C)c1ccccc21",
        "10-hydroxycanthin-6-one": "Oc1ccc2c(c1)c1ccnc3ccc(=O)n2c13",
        "Leurosine": "CC[C@]12C[N@@]3C[C@@H](C[C@@](C(=O)OC)(c4[nH]c5ccccc5c4CC3)c3cc4c(cc3OC)N(C)[C@@H]3[C@]44CCN5CC=C[C@](CC)([C@@H]45)[C@@H](OC(C)=O)[C@]3(O)C(=O)OC)[C@H]1O2",
        "(-)-folicanthine": "[H][C@]12N(C)CC[C@]1(c1ccccc1N2C)[C@@]12CCN(C)[C@]1([H])N(C)c1ccccc21",
        "fumiquinazoline A": "[H][C@]12N[C@@H](C)C(=O)N1c1ccccc1[C@@]2(O)C[C@@H]1C(=O)N[C@@H](C)c2nc3ccccc3c(=O)n12",
        "indole-3-carbaldehyde": "C12=C(NC=C1C([H])=O)C=CC=C2",
        "lycorenan": "[H][C@@]12CC=C3CCN[C@@]3([H])[C@]1([H])c1ccccc1CO2",
        # Some molecules known to be (wrongly) classified:
        "(4-methoxyindol-3-ylmethyl)isothiocyanate": "COC1=CC=CC2=C1C(=CN2)CN=C=S",
        "Gln-Trp-Leu": "O=C1N2[C@H](C(=O)N[C@H](C(=O)N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@@H]1CC3=CC=C(O)C=C3)=O)CC(C)C)=O)O)=O)CC=4C5=C(C=CC=C5)NC4)=O)C(O)C6=CC=C(OC)C=C6)=O)C(O)C7=CC=CC=C7)C[C@H](C2)O",
        "Trp-Val-Pro": "O=C(N[C@@H](CC(C)C)C(O)=O)[C@@H](NC(=O)[C@@H](N)CCC(=O)N)CC=1C=2C(NC1)=CC=CC2"
    }
    
    for name, smi in test_examples.items():
        result, reason = is_indole_alkaloid(smi)
        print(f"{name}: {result}\nReason: {reason}\n")