"""
Classifies: CHEBI:38958 indole alkaloid
"""
#!/usr/bin/env python
"""
Classifies: Indole Alkaloid – An alkaloid containing an indole skeleton.

An indole alkaloid must contain an indole core,
that is, a benzene ring fused with a pyrrole ring. Since many indoles are N–substituted,
we use two SMARTS patterns. In addition, many molecules that are peptides also contain an indole
(e.g., tryptophan residues) – so we add a filter that rejects molecules that contain multiple (2+) amide bonds,
which are far more common in peptides than in small-molecule alkaloids.
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    
    We first parse the SMILES and ensure the molecule is valid. Then a peptide filter is run 
    by counting amide bonds (using the substructure "[NX3][CX3](=O)"). Finally, we check if the molecule 
    matches one of two indole SMARTS queries (one for a typical [nH] indole and one for an N‐substituted indole).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an indole alkaloid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Molecule could not be sanitized"
    
    # Peptide filter: count amide bonds. We consider a substructure "[NX3][CX3](=O)".
    amide_smarts = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 2:
        # Many peptides contain multiple amide bonds.
        return False, "Molecule appears to contain multiple amide bonds (peptidic in nature)"
    
    # Define two SMARTS queries for the indole skeleton.
    # Pattern 1: Classic indole ring with a pyrrole NH.
    indole_smarts1 = Chem.MolFromSmarts("c1cc2c(c1)[nH]c(c2)")
    # Pattern 2: N–substituted indole, which also retains the fused bicyclic structure.
    indole_smarts2 = Chem.MolFromSmarts("c1cc2c(c1)[n]c(c2)")

    # Check for either substructure.
    if mol.HasSubstructMatch(indole_smarts1) or mol.HasSubstructMatch(indole_smarts2):
        return True, ("Molecule contains an indole core (benzene fused to a pyrrole ring); "
                      "thus it can be classified as an indole alkaloid.")
    
    return False, "Molecule does not contain an indole core characteristic of indole alkaloids."

# Example usage:
if __name__ == "__main__":
    # Test some examples:
    test_examples = {
        "staurosporine": "CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)N1C3=C(C=CC=C3)C3=C1C1=C(C4=C(C=CC=C4)N21)C1=C3CNC1=O",
        "indole-3-carbaldehyde": "C12=C(NC=C1C([H])=O)C=CC=C2",
        "Gln-Trp-Leu (peptide)": "O=C(N[C@@H](CC(C)C)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC=2C=3C(NC2)=CC=CC3)C(C)C"
    }
    
    for name, smi in test_examples.items():
        result, reason = is_indole_alkaloid(smi)
        print(f"{name}: {result}\nReason: {reason}\n")