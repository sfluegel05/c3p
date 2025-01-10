"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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
    
    # Define peptide bond pattern (amide bond between amino acids)
    peptide_bond = Chem.MolFromSmarts("N[C;H1,C](=O)")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond)
    
    if not peptide_bonds:
        return False, "No peptide bonds found"
    
    # Identify amino acid residues by looking for alpha carbon patterns
    alpha_carbon = Chem.MolFromSmarts("[C;H1,H2]([C;H2,H3])[C](=O)N")
    residues = mol.GetSubstructMatches(alpha_carbon)
    num_residues = len(residues)
    
    if num_residues < 2:
        return False, f"Only {num_residues} amino acid residue found, need at least 2"
    elif num_residues > 20:
        return False, f"{num_residues} amino acid residues found, exceeds typical oligopeptide length"
    
    # Check if all peptide bonds are connected to form a peptide chain
    # Create a graph to check connectivity
    peptide_bond_indices = [bond[1] for bond in peptide_bonds]  # Get bonded carbonyl carbons
    residue_atoms = [res[0] for res in residues]  # Get alpha carbons
    # Check if peptide bonds connect the residues sequentially
    connected = True
    for i in range(len(residue_atoms)-1):
        path = Chem.GetShortestPath(mol, residue_atoms[i], residue_atoms[i+1])
        if not any(atom_idx in peptide_bond_indices for atom_idx in path):
            connected = False
            break
    if not connected:
        return False, "Amino acids are not connected via peptide bonds"
    
    return True, f"Oligopeptide with {num_residues} amino acid residues"