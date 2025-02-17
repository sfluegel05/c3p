"""
Classifies: CHEBI:38958 indole alkaloid
"""
#!/usr/bin/env python
"""
Classifies: Indole Alkaloid – An alkaloid containing an indole skeleton.
An indole alkaloid must contain an indole core (a benzene ring fused to a pyrrole ring). 
This program checks for the presence of an indole moiety using two SMARTS patterns – 
one for an indole with a free NH and one for an N‐substituted indole.
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    
    The function parses the SMILES string, sanitizes the molecule, and then checks if the
    molecule contains an indole core. An indole core is defined as a benzene ring fused to 
    a pyrrole ring. Two SMARTS patterns are used: one for the classic indole (free NH) and one 
    for an N–substituted indole.
    
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
    
    # Sanitize the molecule.
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return False, "Molecule could not be sanitized"
    
    # We previously attempted to filter out peptides by counting amide bonds.
    # However, many indole alkaloids contain amide groups, so we no longer apply that filter.
    
    # Define SMARTS patterns for the indole core.
    # Pattern 1: A classic indole ring (benzene ring fused with a pyrrole having a free NH).
    indole_smarts1 = Chem.MolFromSmarts("c1ccc2[nH]c(c2)c1")
    # Pattern 2: N-substituted indole (fused structure with the nitrogen substituted).
    indole_smarts2 = Chem.MolFromSmarts("c1ccc2c(c1)[n]cc2")
    
    # Check for a match to either SMARTS pattern.
    if mol.HasSubstructMatch(indole_smarts1) or mol.HasSubstructMatch(indole_smarts2):
        return True, ("Molecule contains an indole core (benzene fused to a pyrrole ring) "
                      "and can be classified as an indole alkaloid.")
    
    return False, "Molecule does not contain an indole core characteristic of indole alkaloids."

# Example usage when running as a script.
if __name__ == "__main__":
    # Some indole alkaloid examples (the SMILES strings given in the prompt).
    test_examples = {
        "staurosporine": "CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)N1C3=C(C=CC=C3)C3=C1C1=C(C4=C(C=CC=C4)N21)C1=C3CNC1=O",
        "3alpha(S)-strictosidine": "[H][C@@]1(C[C@H]2[C@@H](C=C)[C@H](O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)OC=C2C(=O)OC)NCCC2=C1NC1=C2C=CC=C1",
        "folicanthine": "CN1CCC2(C1N(C)c1ccccc21)C12CCN(C)C1N(C)c1ccccc21",
        "10-hydroxycanthin-6-one": "Oc1ccc2c(c1)c1ccnc3ccc(=O)n2c13",
        "Leurosine": "CC[C@]12C[N@@]3C[C@@H](C[C@@](C(=O)OC)(c4[nH]c5ccccc5c4CC3)c3cc4c(cc3OC)N(C)[C@@H]3[C@]44CCN5CC=C[C@](CC)([C@@H]45)[C@@H](OC(C)=O)[C@]3(O)C(=O)OC)[C@H]1O2",
        "(-)-folicanthine": "[H][C@]12N(C)CC[C@]1(c1ccccc1N2C)[C@@]12CCN(C)[C@]1([H])N(C)c1ccccc21",
        "fumiquinazoline A": "[H][C@]12N[C@@H](C)C(=O)N1c1ccccc1[C@@]2(O)C[C@@H]1C(=O)N[C@@H](C)c2nc3ccccc3c(=O)n12",
        "indole-3-carbaldehyde": "C12=C(NC=C1C([H])=O)C=CC=C2",
        "lycorenan": "[H][C@@]12CC=C3CCN[C@@]3([H])[C@]1([H])c1ccccc1CO2"
    }
    
    # Test the classifier for each sample.
    for name, smi in test_examples.items():
        result, reason = is_indole_alkaloid(smi)
        print(f"{name}: {result}\nReason: {reason}\n")