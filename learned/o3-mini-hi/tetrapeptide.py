"""
Classifies: CHEBI:48030 tetrapeptide
"""
#!/usr/bin/env python3
"""
Classifies: Tetrapeptide
Definition: Any molecule that contains four amino‐acid residues connected by peptide linkages.
A tetrapeptide (linear) should possess three backbone peptide bonds.
This function uses a heuristic: it looks for the substructure "C(=O)N" and then filters
the matches to ensure that the nitrogen is attached to a chiral carbon (i.e. the alpha carbon).
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is defined as any molecule that contains four amino acid residues connected by peptide linkages.
    For a linear tetrapeptide, there should be three peptide bonds in the main chain.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is (likely) a tetrapeptide, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a generic amide bond: C(=O)N.
    # Note: This substructure might also capture side-chain amides, so we filter further.
    peptide_bond_smarts = "C(=O)N"
    pattern = Chem.MolFromSmarts(peptide_bond_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern used for peptide bond detection"

    matches = mol.GetSubstructMatches(pattern)
    # Use a set to keep unique peptide bond candidates (by their atom indices).
    qualifying_bonds = set()

    # For each match, check if the N atom is attached to a chiral (alpha carbon) atom.
    # In peptides many residues (except glycine) are defined with a chiral center ([C@@H] or [C@H])
    # which is usually the alpha carbon.
    for match in matches:
        # match returns a tuple of atom indices corresponding to the pattern atoms (C, then N)
        c_idx, n_idx = match
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Look for a neighbor of the N atom which is a carbon with an explicit chiral tag.
        has_chiral_neighbor = False
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetIdx() == c_idx:
                continue  # Skip the carbonyl carbon already matched.
            if neighbor.GetSymbol() == "C" and neighbor.GetChiralTag() != rdchem.ChiralType.CHI_UNSPECIFIED:
                has_chiral_neighbor = True
                break
        if has_chiral_neighbor:
            qualifying_bonds.add((c_idx, n_idx))
    
    count = len(qualifying_bonds)
    
    # A linear tetrapeptide should have three backbone peptide bonds.
    if count == 3:
        return True, "Found 3 backbone peptide bonds consistent with a tetrapeptide structure."
    else:
        return False, f"Found {count} backbone peptide bonds; expected 3 for a tetrapeptide."

# Example usage (uncomment to test):
# smiles_example = "C[C@@H](O)[C@H](NC(=O)[C@@H](N)CCCCN)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O"  # Lys-Thr-Trp-Tyr
# result, reason = is_tetrapeptide(smiles_example)
# print(result, reason)