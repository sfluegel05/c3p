"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: Polypeptide, defined as a peptide containing ten or more amino acid residues.
A polypeptide is here heuristically identified by counting peptide bonds.
For a linear peptide, residue count = (number of amide bonds) + 1
For a cyclic peptide (where terminal groups are connected), residue count = (number of amide bonds).
This implementation uses simple RDKit substructure searches to count the characteristic “C(=O)N”
pattern and to judge whether the molecule is linear (with free N‑ and C‑termini) or cyclic.
Note: This is a heuristic approach that may not capture all edge‐cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (≥10 amino acid residues) based on its SMILES string.
    Uses a heuristic by counting the number of peptide bonds (C(=O)N).
    For a linear peptide, number of residues = number of peptide bonds + 1;
    for a cyclic peptide, the number of residues is taken as equal to the number of peptide bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a polypeptide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for an amide bond: C(=O)N
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    n_amide = len(amide_matches)
    
    # Heuristic: Check for free termini which are typical for linear peptides.
    # A free N-terminal amine (primary amine) is roughly matched by [NH2]
    # and a free C-terminal acid by C(=O)[OH].
    free_amine_pattern = Chem.MolFromSmarts("[NH2]")
    free_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    has_free_amine = mol.HasSubstructMatch(free_amine_pattern)
    has_free_acid = mol.HasSubstructMatch(free_acid_pattern)
    
    # If both free amine and free acid are found, treat the peptide as linear.
    if has_free_amine and has_free_acid:
        residue_count = n_amide + 1
        peptide_type = "linear"
    else:
        # Otherwise, assume a cyclic (or “blocked‐termini”) peptide.
        residue_count = n_amide  # In cyclic peptides, the number of peptide bonds equals the residue count.
        peptide_type = "cyclic or blocked termini"
    
    # Classify as polypeptide if residue_count is 10 or more.
    if residue_count >= 10:
        return True, f"Detected {residue_count} amino acid residues (peptide appears {peptide_type})."
    else:
        return False, f"Detected only {residue_count} amino acid residues (need at least 10)."
        
# Example usage:
if __name__ == "__main__":
    # Some test SMILES strings (these are among the examples provided)
    smiles_list = [
        "CC[C@H](C)[C@H](NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCNC(=O)c1ccccc1C(F)(F)F)NC(=O)[C@H](CC(O)=O)NC(C)=O)[C@@H](C)O)C(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(O)=O",
        "CC[C@H](C)[C@H](NC(=O)[C@H](CO)NC(=O)[C@H](CC1=CC=C(O)C=C1)NC(=O)[C@H](CO)NC(C)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(N)=O)C(=O)NCC(=O)NCC(=O)N1CCC[C@H]1C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](C)C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(=O)N[C@@H](CO)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(N)=O"
    ]
    for sm in smiles_list:
        result, reason = is_polypeptide(sm)
        print("SMILES:", sm)
        print("Is polypeptide?", result)
        print("Reason:", reason, "\n")