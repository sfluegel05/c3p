"""
Classifies: CHEBI:60834 beta-carbolines
"""
#!/usr/bin/env python
"""
Classifies: beta-carbolines (any pyridoindole containing a beta-carboline skeleton and their hydrogenated derivatives)
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string) is a beta-carboline or a hydrogenated derivative.
    The definition used here is that of any pyridoindole containing a beta-carboline skeleton.
    This function uses substructure matching with a SMARTS pattern intended to capture the typical beta-carboline tricyclic core.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a beta-carboline skeleton, False otherwise.
        str: Reason for classification.
    """
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for the beta-carboline (pyrido[3,4-b]indole) skeleton.
    # This pattern looks for a tricyclic system where an indole (benzene fused to a pyrrole with a [nH]) is fused to a pyridine-type ring.
    # Note: This pattern is formulated in aromatic form and will detect many beta-carboline structures,
    # though highly hydrogenated derivatives may sometimes escape detection.
    beta_carboline_smarts = "c1ccc2c(c1)[nH]c3c(n2)cccc3"
    query = Chem.MolFromSmarts(beta_carboline_smarts)
    
    # Check for substructure match
    if mol.HasSubstructMatch(query):
        return True, "Molecule contains the beta-carboline (pyridoindole) skeleton."
    else:
        # In some cases (especially hydrogenated derivatives) the aromatic flags might not be set.
        # One can try to kekulize the molecule (which often reassigns aromaticity) and check again.
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception:
            pass
        if mol.HasSubstructMatch(query):
            return True, "Molecule contains the beta-carboline skeleton (detected after kekulization)."
    
    return False, "Molecule does not appear to contain a beta-carboline (pyridoindole) skeleton."

# Example usage:
if __name__ == "__main__":
    # Test using harmine (a well-known beta-carboline) SMILES:
    harmine_smiles = "COc1ccc2c(c1)[nH]c1c(C)nccc21"
    result, reason = is_beta_carbolines(harmine_smiles)
    print("Harmine:", result, reason)