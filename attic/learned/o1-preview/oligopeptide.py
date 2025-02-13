"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: oligopeptide
"""
from rdkit import Chem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing a relatively small number of amino acids (typically 2 to 20 residues).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all amide bonds not in rings (peptide bonds)
    amide_bond_smarts = "[NX3][CX3](=O)"
    amide_bond_pattern = Chem.MolFromSmarts(amide_bond_smarts)
    amide_matches = mol.GetSubstructMatches(amide_bond_pattern)
    
    num_peptide_bonds = 0
    for match in amide_matches:
        n_idx = match[0]
        c_idx = match[1]
        # Check if the bond between N and C is not in a ring
        bond = mol.GetBondBetweenAtoms(n_idx, c_idx)
        if bond and not bond.IsInRing():
            num_peptide_bonds += 1
    
    if num_peptide_bonds == 0:
        return False, "No peptide bonds found"
    
    # Estimate number of amino acid residues
    num_residues = num_peptide_bonds + 1  # For linear peptides

    if num_residues < 2:
        return False, f"Only {num_residues} amino acid residue found, need at least 2"
    elif num_residues > 20:
        return False, f"{num_residues} amino acid residues found, exceeds typical oligopeptide length"
    
    return True, f"Oligopeptide with approximately {num_residues} amino acid residues"