"""
Classifies: CHEBI:48030 tetrapeptide
"""
Below is one idea. I first analyzed the previous code’s problems. It was based solely on a generic “C(=O)N” pattern and then a chiral‐neighbor check. In practice many compounds (for example, oligosaccharides with acetamido groups or other amide bonds) can also have that pattern, and some “true” tetrapeptides (especially when protected or modified) can have extra amide bonds. One strategy is to “zoom in” on the full peptide‐backbone motif: namely, look for a peptide bond together with an immediately attached (chiral) α‐carbon. In a “linear” tetrapeptide there should be three backbone linkages (which connect four α–carbons). In addition, one can check that the overall molecular weight and overall nitrogen count fall into a reasonable range. (Tetrapeptides tend to be—not very large molecules, and they usually have one amino‐nitrogen per residue; protecting groups and extra functionalities may change these numbers but we can at least “penalize” very heavy molecules.) 

The program below uses two SMARTS patterns (one for each possible explicit chirality of the α–carbon) in order to “capture” amide bonds that are immediately followed by an α–carbon (i.e. it looks for a fragment like “C(=O)N[C@H]” or “C(=O)N[C@@H]”). For every match we record the “peptide bond” (as a pairing of the carbonyl carbon with the amide nitrogen) and then require that exactly three such bonds are found. Finally, we check that the overall molecular weight is in a typical tetrapeptide range and that the total number of nitrogen atoms is reasonable. (These additional filters do not “prove” that a molecule is a tetrapeptide but may reduce some of the many false hits.) 

Below is one possible solution in a single Python code block:

------------------------------------------------

#!/usr/bin/env python3
"""
Classifies: Tetrapeptide
Definition: Any molecule that contains four amino‐acid residues connected by peptide linkages.
For a linear tetrapeptide, exactly three backbone peptide bonds (linking four α–carbons) are expected.
The heuristic implemented here first searches for peptide‐bond motifs
that include an attached, explicitly stereochemically–defined α–carbon (i.e. "C(=O)N[C@H]" or "C(=O)N[C@@H]").
Then additional filters (molecular weight and nitrogen count) are applied.
Note: This heuristic is not perfect and may mis‐classify some molecules.
"""
from rdkit import Chem
from rdkit.Chem import rdchem, rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A linear tetrapeptide (or one that is “backbone–tetrapeptidic”) should have three peptide bonds that 
    join four amino–acid residues. Here we look for the chain fragment C(=O)N[C@H] (or with C@@H),
    then add extra filters on molecular weight and overall nitrogen count.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is (likely) a tetrapeptide, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for peptide bonds with an attached (chiral) α–carbon.
    # These 4–atom patterns capture the fragment: [C:carbonyl](=O)[N][C@H] or [C:carbonyl](=O)[N][C@@H]
    pat1 = Chem.MolFromSmarts("C(=O)N[C@H]")
    pat2 = Chem.MolFromSmarts("C(=O)N[C@@H]")
    if pat1 is None or pat2 is None:
        return False, "Error creating SMARTS patterns"
    
    peptide_bonds = set()  # will hold (carbonyl carbon index, amide nitrogen index) pairs
    # Find matches with the first pattern.
    for match in mol.GetSubstructMatches(pat1):
        # Expected match length is 4: (carbonyl C, carbonyl O, amide N, α–carbon)
        if len(match) != 4: 
            continue
        c_idx = match[0]   # carbonyl carbon
        n_idx = match[2]   # amide nitrogen
        peptide_bonds.add((c_idx, n_idx))
        
    # Do the same for the second pattern.
    for match in mol.GetSubstructMatches(pat2):
        if len(match) != 4:
            continue
        c_idx = match[0]
        n_idx = match[2]
        peptide_bonds.add((c_idx, n_idx))
    
    n_bonds = len(peptide_bonds)
    
    # For a linear tetrapeptide we expect exactly three backbone peptide bonds.
    if n_bonds != 3:
        return False, f"Found {n_bonds} peptide bond(s) matching the motif; expected 3 for a linear tetrapeptide."
    
    # Additional filters: check that the molecular weight is in a typical range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (250 <= mol_wt <= 1500):
        return False, f"Molecular weight {mol_wt:.1f} Da out of typical tetrapeptide range (250-1500 Da)."
    
    # And check that the number of nitrogen atoms is reasonable.
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    # Typical tetrapeptides have at least one backbone N per residue but usually not far above 10 (unless extra N–rich protecting groups)
    if n_nitrogen < 4 or n_nitrogen > 12:
        return False, f"Found {n_nitrogen} nitrogen atom(s), which is not typical for a tetrapeptide."
    
    # If passed all checks we tentatively classify as tetrapeptide.
    return True, "Found 3 backbone peptide bonds with adjacent α–carbons and overall properties consistent with a tetrapeptide."

# Uncomment the following lines to run a test (example of Lys-Thr-Trp-Tyr tetrapeptide):
# test_smiles = "C[C@@H](O)[C@H](NC(=O)[C@@H](N)CCCCN)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O"
# result, reason = is_tetrapeptide(test_smiles)
# print(result, reason)

# Note: This heuristic may miss or wrongly classify some molecules (especially those having extra amide bonds
# from protecting groups or modifications, or macrocyclic peptides). Further refinements could include
# a graph search for a contiguous chain of four α–carbons connected by peptide bonds.
    
if __name__ == '__main__':
    # Testing on a sample (the Lys-Thr-Trp-Tyr tetrapeptide)
    test_smiles = "C[C@@H](O)[C@H](NC(=O)[C@@H](N)CCCCN)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O"
    result, reason = is_tetrapeptide(test_smiles)
    print(result, reason)

----------------------------------------

In this solution the function first looks for “peptide‐bond fragments” that include an immediately attached chiral α–carbon. Then it checks that exactly three such bonds were found (a linear chain of four residues). Finally, molecular‐weight and nitrogen–count filters are used. This approach (like the previous one) is heuristic – it may mis‐classify some modified, cyclic, or heavily protected species, but it should improve on the many false classifications seen before.