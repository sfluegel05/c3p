"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a short chain of amino acids linked by peptide bonds.

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

    # Look for peptide bonds (-C(=O)N-)
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_matches)
    if num_peptide_bonds < 1:
        return False, "Found no peptide bonds"
    
    amino_acid_residues = 0
    for match in peptide_matches:
        c_atom_index = match[0]
        n_atom_index = match[1]
        c_atom = mol.GetAtomWithIdx(c_atom_index)
        n_atom = mol.GetAtomWithIdx(n_atom_index)

        # Find alpha carbon
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 4:
                alpha_carbon = neighbor
                for alpha_neighbor in alpha_carbon.GetNeighbors():
                     if alpha_neighbor.GetIdx() == n_atom_index:
                         amino_acid_residues += 1
                         break;
    
    if amino_acid_residues < 2:
        return False, f"Found {amino_acid_residues} amino acid residues, at least 2 required"
    
    # Check for number of residues
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 4 or num_atoms > 500: #Rough sanity check
        return False, f"Molecule has {num_atoms}, which is outside of the reasonable size for oligopeptides"
    
    return True, "Has peptide bonds and amino acid residues, within size range of a typical oligopeptide"