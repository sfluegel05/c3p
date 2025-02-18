"""
Classifies: CHEBI:60834 beta-carbolines
"""
#!/usr/bin/env python
"""
Classifies: beta-carbolines 
Definition: Any pyridoindole containing a beta-carboline skeleton and their hydrogenated derivatives.
Note: Beta-carbolines are built on a fused indole-pyridine (pyridoindole) system. The aromatic form is typical 
(e.g. harmine) while many of the “dihydro” or hydrogenated variants have less aromaticity. Because substructure 
searches are approximations, we use two SMARTS queries:
   • The first (arom_query) looks for the fully aromatic beta-carboline core.
   • The second (hydro_query) is intended to catch some hydrogenated (dihydro) beta-carboline variants.
This approach is an approximation and may still lead to false positives or negatives.
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string) is a beta-carboline or its hydrogenated derivative.
    It uses two SMARTS queries:
      1. arom_query: Looks for an aromatic beta-carboline core, defined as a pyrido[3,4-b]indole ring system.
         (SMARTS: "c1ccc2c(c1)[nH]c3ccccc23")
      2. hydro_query: Looks for a partially/hydrogenated variant.
         Here we use a pattern that roughly captures a fused ring system having an indole-like ring (with [nH])
         fused to a saturated or partly saturated six-membered ring.
         (SMARTS: "C1CN2CCc3c[nH]c(c3)C2C1")
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple of classification (True if beta-carboline core is found) and a reason string.
    """
    # Parse SMILES string into molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # --- First SMARTS: aromatic beta-carboline ---
    # This pattern looks for a pyrido[3,4-b]indole (the indole ring fused to a pyridine ring):
    arom_smarts = "c1ccc2c(c1)[nH]c3ccccc23"
    arom_query = Chem.MolFromSmarts(arom_smarts)
    if arom_query is None:
        return False, "Error in constructing aromatic SMARTS."
    
    if mol.HasSubstructMatch(arom_query):
        return True, "Molecule contains an aromatic beta-carboline (pyridoindole) skeleton."
    
    # Sometimes aromatic perception (or bond order assignment) may cause the SMARTS match to fail.
    # Try kekulizing and re-checking.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass
    if mol.HasSubstructMatch(arom_query):
        return True, "Molecule contains an aromatic beta-carboline skeleton (detected after kekulization)."
    
    # --- Second SMARTS: hydrogenated (dihydro) beta-carboline ---
    # This query is less strict about aromaticity. It attempts to catch a fused system where an indole-like ring
    # (with an NH) is fused to a second ring that is at least partially saturated.
    # Note: This pattern will not capture all possible variants and may mis‐match cases.
    hydro_smarts = "C1CN2CCc3c[nH]c(c3)C2C1"
    hydro_query = Chem.MolFromSmarts(hydro_smarts)
    if hydro_query is None:
        return False, "Error in constructing hydrogenated SMARTS."
    
    if mol.HasSubstructMatch(hydro_query):
        return True, "Molecule contains a hydrogenated (dihydro) beta-carboline skeleton."
    
    # If neither query finds a match, then we do not classify this as a beta-carboline.
    return False, "Molecule does not appear to contain a beta-carboline (pyridoindole) skeleton."

# Example usage (for testing) when running the script directly.
if __name__ == "__main__":
    test_examples = {
        "harmine (aromatic beta-carboline)": "COc1ccc2c(c1)[nH]c3ccccc23",
        "2-Methyl-6-methoxy-1,2,3,4-tetrahydro-beta-carboline": "O(C1=CC=2C3=C(NC2C=C1)CN(CC3)C)C",
    }
    for name, smi in test_examples.items():
        result, reason = is_beta_carbolines(smi)
        print(f"{name}: {result} ({reason})")