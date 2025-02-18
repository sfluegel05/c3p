"""
Classifies: CHEBI:60834 beta-carbolines
"""
#!/usr/bin/env python
"""
Classifies: beta-carbolines 
Definition: Any pyridoindole containing a beta-carboline skeleton and their hydrogenated derivatives.

This method uses two SMARTS queries:
  • One query (arom_query) is geared toward fully aromatic beta‐carbolines. It requires, for example,
    the presence of an indole fragment having an [nH] atom fused to a pyridine‐type ring.
  • The second query (hydro_query) is intended to catch partially/saturated (“dihydro”) beta‐carbolines.
Because “rule‐based” substructure searches can sometimes over‐ or under‐match,
we also try a kekulization attempt if the first search fails.
Note: This set of patterns is an attempt to improve on previous failures and may still not perfectly separate all edge cases.
"""

from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string) is considered a beta-carboline (or a hydrogenated derivative).
    The definition is: any pyridoindole containing a beta-carboline skeleton.
    
    This function uses two SMARTS queries:
      - arom_query: looks for a fully aromatic beta-carboline core.
      - hydro_query: looks for a likely hydrogenated (dihydro) beta-carboline variant.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple of a boolean (True if beta-carboline core is found, False otherwise)
                     and a reason string.
    """
    # Parse the input molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS pattern for fully aromatic beta-carboline.
    # This pattern requires a fully aromatic fused system where the indole part has an [nH] and is fused to another ring.
    arom_smarts = "c1ccc2c(c1)[nH]c3c2cccc3"
    arom_query = Chem.MolFromSmarts(arom_smarts)
    if arom_query is None:
        return False, "Error in constructing aromatic SMARTS."
    
    # Check for the aromatic beta-carboline skeleton
    if mol.HasSubstructMatch(arom_query):
        return True, "Molecule contains an aromatic beta-carboline (pyridoindole) skeleton."
    
    # In case aromaticity assignment (or the bond orders) did not favor detection,
    # try to kekulize (which reassigns bond orders) and check again.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass  # if kekulization fails, continue
    
    if mol.HasSubstructMatch(arom_query):
        return True, "Molecule contains an aromatic beta-carboline skeleton (detected after kekulization)."
    
    # Define a second SMARTS to catch hydrogenated (e.g. dihydro or tetrahydro) beta-carboline variants.
    # This pattern is less strict about aromaticity. It attempts to match a fused system in which
    # a six-membered ring with one nitrogen is fused to an indole-like motif.
    # (This pattern is an approximation and may not be perfect.)
    hydro_smarts = "C1CC2=CC=CC=C2NC1"
    hydro_query = Chem.MolFromSmarts(hydro_smarts)
    if hydro_query is None:
        return False, "Error in constructing hydrogenated SMARTS."
    
    if mol.HasSubstructMatch(hydro_query):
        return True, "Molecule contains a hydrogenated beta-carboline skeleton."
    
    # If neither aromatic nor hydrogenated queries match, then do not classify as beta-carboline.
    return False, "Molecule does not appear to contain a beta-carboline (pyridoindole) skeleton."

# Example usage (for testing)
if __name__ == "__main__":
    # Test examples: harmine (aromatic beta-carboline) and one dihydro-type.
    examples = {
        "harmine": "COc1ccc2c(c1)[nH]c1c(C)nccc21",
        "2-Methyl-6-methoxy-1,2,3,4-tetrahydro-beta-carboline": "O(C1=CC=2C3=C(NC2C=C1)CN(CC3)C)C"
    }
    for name, smi in examples.items():
        result, reason = is_beta_carbolines(smi)
        print(f"{name}: {result} ({reason})")