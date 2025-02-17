"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: Oligopeptide – a peptide containing a relatively small number of amino acids.
This heuristic function checks for amide bonds (the peptide bonds) and estimates the residue count.
We assume that an oligopeptide has between 2 and 10 amino acids (i.e., 1 to 9 amide bonds).
Additionally, the presence of an alpha-carbon unit [C@H](N)C(=O) or [C@@H](N)C(=O) is verified.
Note: This heuristic will not catch every nuance.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide from its SMILES string.
    A peptide is characterized by amide bonds (–C(=O)N–), and we estimate the number of residues 
    (number of peptide bonds + 1) should be relatively small (here assumed <=10).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an oligopeptide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for an amide (peptide) bond.
    peptide_bond_smarts = Chem.MolFromSmarts("C(=O)N")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_smarts)
    n_peptide_bonds = len(peptide_bonds)
    
    # Heuristic: if there are no peptide bonds, it's not a peptide.
    if n_peptide_bonds == 0:
        return False, "No peptide (amide) bonds found"
    
    # Estimate number of amino acid residues (peptide bonds + 1).
    n_residues = n_peptide_bonds + 1
    if n_residues > 10:
        return False, f"Found {n_residues} amino acid residues which is too many for an oligopeptide"
    
    # Check for at least one typical amino acid unit with an alpha carbon using chirality tags.
    # We will look for patterns of the form [C@H](N)C(=O) and [C@@H](N)C(=O)
    aa_pattern1 = Chem.MolFromSmarts("[C@H](N)C(=O)")
    aa_pattern2 = Chem.MolFromSmarts("[C@@H](N)C(=O)")
    match1 = mol.HasSubstructMatch(aa_pattern1)
    match2 = mol.HasSubstructMatch(aa_pattern2)
    if not (match1 or match2):
        # It might still be a peptide if the terminal unit lost chirality, but likely a peptide backbone is absent.
        return False, "No typical alpha-amino acid unit (with chiral center) found"
    
    # Optionally we can calculate some additional properties (rotatable bonds, molecular weight)
    # that are common in peptides.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Give some extra info in the reason
    reason = (f"Detected {n_peptide_bonds} peptide bonds (≈{n_residues} residues), "
              f"{n_rotatable} rotatable bonds, and MW of {mol_wt:.1f} Da. "
              "This is consistent with an oligopeptide.")
    
    return True, reason