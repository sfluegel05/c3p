"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: CHEBI:33709 amino acid
A carboxylic acid containing one or more amino groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for amino groups (primary, secondary, tertiary)
    amino_pattern = Chem.MolFromSmarts("[N;!H0;!H1;!H2;!H3;!$(NC=O)]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if not amino_matches:
        return False, "No amino group found"
    
    # Exclude peptides and other molecules with multiple amide bonds
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) > 1:
        return False, "Contains multiple amide bonds, likely a peptide"
    
    # Check for single chiral center (typical of amino acids)
    chiral_centers = AllChem.FindMolChiralCenters(mol, includeUnspecified=True)
    if len(chiral_centers) != 1:
        return False, "Does not have a single chiral center"
    
    # Check for specific arrangement of carboxylic acid and amino groups
    chiral_center = chiral_centers[0]
    chiral_atom = mol.GetAtomWithIdx(chiral_center)
    neighbors = [mol.GetAtomWithIdx(n) for n in chiral_atom.GetNeighbors()]
    
    carboxyl_neighbor = None
    amino_neighbor = None
    for neighbor in neighbors:
        if neighbor.GetSymbol() == "O" and neighbor.GetFormalCharge() == -1:
            carboxyl_neighbor = neighbor
        elif neighbor.GetSymbol() == "N" and neighbor.GetFormalCharge() == 0:
            amino_neighbor = neighbor
    
    if carboxyl_neighbor is None or amino_neighbor is None:
        return False, "Incorrect arrangement of carboxylic acid and amino groups"
    
    # Handle potential zwitterions by considering alternative resonance structures
    zwitterion_pattern = Chem.MolFromSmarts("[N+;H2,H1;!$(NC=O)]")
    if mol.HasSubstructMatch(zwitterion_pattern):
        zwitterion_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))
        amino_matches = zwitterion_mol.GetSubstructMatches(amino_pattern)
        if amino_matches:
            return True, "Amino acid (zwitterionic form)"
    
    return True, "Amino acid"