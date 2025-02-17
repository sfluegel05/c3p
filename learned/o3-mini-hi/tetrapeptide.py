"""
Classifies: CHEBI:48030 tetrapeptide
"""
#!/usr/bin/env python3
"""
Classifies: Tetrapeptide
Definition: Any molecule that contains four amino‐acid residues connected by peptide linkages.
For a linear tetrapeptide, there should be three backbone peptide bonds.
This function uses a heuristic: it looks for the substructure "C(=O)N" (which is actually 3 atoms)
and then filters the matches to ensure that the nitrogen is attached to a chiral (alpha) carbon.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is defined as any molecule that contains four amino acid residues connected by peptide linkages.
    For a linear tetrapeptide, there should be exactly three backbone peptide bonds.
    
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
    
    # Define a SMARTS pattern for a generic amide (peptide) bond.
    # Note that "C(=O)N" will match three atoms: C, O, and N.
    peptide_bond_smarts = "C(=O)N"
    pattern = Chem.MolFromSmarts(peptide_bond_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern used for peptide bond detection"
    
    # Find all substructure matches for the peptide bond.
    matches = mol.GetSubstructMatches(pattern)
    
    qualifying_bonds = set()
    
    # For each substructure match, ensure that the nitrogen is attached to a chiral (alpha) carbon.
    for match in matches:
        # Ensure that exactly three atoms were matched (C, O, N).
        if len(match) != 3:
            continue  # skip unexpected matches
        c_idx, o_idx, n_idx = match
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check if N has a neighbor (other than the carbonyl carbon) which is a carbon with an explicit chiral tag.
        has_chiral_neighbor = False
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetIdx() == c_idx:
                continue  # skip the already matched carbonyl carbon
            if neighbor.GetSymbol() == "C" and neighbor.GetChiralTag() != rdchem.ChiralType.CHI_UNSPECIFIED:
                has_chiral_neighbor = True
                break
        if has_chiral_neighbor:
            qualifying_bonds.add((c_idx, n_idx))
    
    count = len(qualifying_bonds)
    
    # A linear tetrapeptide should have exactly three backbone peptide bonds.
    if count == 3:
        return True, "Found 3 backbone peptide bonds consistent with a tetrapeptide structure."
    else:
        return False, f"Found {count} backbone peptide bonds; expected 3 for a tetrapeptide."

# Example usage (uncomment the following lines to test the function):
# smiles_example = "C[C@@H](O)[C@H](NC(=O)[C@@H](N)CCCCN)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O"  # Lys-Thr-Trp-Tyr tetrapeptide example
# result, reason = is_tetrapeptide(smiles_example)
# print(result, reason)