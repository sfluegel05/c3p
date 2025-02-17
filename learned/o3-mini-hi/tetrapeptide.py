"""
Classifies: CHEBI:48030 tetrapeptide
"""
#!/usr/bin/env python3
"""
Classifies: Tetrapeptide
Definition: Any molecule that contains four amino-acid residues connected by peptide linkages.
For a linear tetrapeptide, exactly three peptide bonds (linking four alpha-carbons) are expected.
This heuristic searches for the fragment "C(=O)N[C@H]" (or with [C@@H]) to count backbone peptide bonds,
and then applies additional filters on molecular weight and nitrogen count.
Note: This heuristic may miss or misclassify molecules with modifications, cyclic peptides, or extra amide bonds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A linear tetrapeptide should have three peptide bonds connecting four alpha-carbons.
    The implementation searches for backbone peptide bond motifs with an explicitly chiral alpha-carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a tetrapeptide, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse the SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns to capture the fragments "C(=O)N[C@H]" and "C(=O)N[C@@H]".
    pat1 = Chem.MolFromSmarts("C(=O)N[C@H]")
    pat2 = Chem.MolFromSmarts("C(=O)N[C@@H]")
    if pat1 is None or pat2 is None:
        return False, "Error creating SMARTS patterns"
    
    peptide_bonds = set()  # Use a set to store unique (carbonyl C, amide N) pairs
    
    # Search for matches using the first pattern.
    for match in mol.GetSubstructMatches(pat1):
        # Expecting four atoms in the fragment: carbonyl C, carbonyl O, amide N, and alpha-C.
        if len(match) != 4:
            continue
        c_idx = match[0]  # carbonyl carbon
        n_idx = match[2]  # amide nitrogen
        peptide_bonds.add((c_idx, n_idx))
    
    # Search for matches using the second pattern.
    for match in mol.GetSubstructMatches(pat2):
        if len(match) != 4:
            continue
        c_idx = match[0]
        n_idx = match[2]
        peptide_bonds.add((c_idx, n_idx))
    
    n_bonds = len(peptide_bonds)
    
    # In a linear tetrapeptide, we expect exactly three peptide bonds.
    if n_bonds != 3:
        return False, "Found {} peptide bond(s) matching the motif; expected 3 for a linear tetrapeptide.".format(n_bonds)
    
    # Check additional properties:
    # Typical tetrapeptides usually lie within a certain molecular weight range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (250 <= mol_wt <= 1500):
        return False, "Molecular weight {:.1f} Da out of typical tetrapeptide range (250-1500 Da).".format(mol_wt)
    
    # Check the number of nitrogen atoms.
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_nitrogen < 4 or n_nitrogen > 12:
        return False, "Found {} nitrogen atom(s), which is not typical for a tetrapeptide.".format(n_nitrogen)
    
    return True, "Found 3 peptide bonds with adjacent chiral alpha-carbons and overall properties consistent with a tetrapeptide."

# For quick testing, uncomment a test case below:
# test_smiles = "C[C@@H](O)[C@H](NC(=O)[C@@H](N)CCCCN)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O"
# result, reason = is_tetrapeptide(test_smiles)
# print(result, reason)

if __name__ == '__main__':
    # Example: Testing using one of the given tetrapeptide SMILES (Lys-Thr-Trp-Tyr)
    test_smiles = "C[C@@H](O)[C@H](NC(=O)[C@@H](N)CCCCN)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O"
    result, reason = is_tetrapeptide(test_smiles)
    print(result, reason)