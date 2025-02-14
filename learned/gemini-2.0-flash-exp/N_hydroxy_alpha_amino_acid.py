"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino acid based on its SMILES string.
    An N-hydroxy-alpha-amino acid is an amino acid in which at least one hydrogen 
    attached to the amino group is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Pattern for alpha carbon (carbon attached to carboxyl group and another carbon)
    alpha_carbon_pattern = Chem.MolFromSmarts("[CH](C(=O)O)[CX4]")
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    if not alpha_carbon_matches:
        return False, "No alpha carbon found"

    # Pattern for N-hydroxy group
    n_hydroxy_pattern_1 = Chem.MolFromSmarts("[NX3][OH1]")  # -NH-OH or -N-OH
    n_hydroxy_pattern_2 = Chem.MolFromSmarts("[NX3]([OH1])[OH1]") # -N(OH)2
    
    found_n_hydroxy = False

    for match in alpha_carbon_matches:
        alpha_carbon = mol.GetAtomWithIdx(match[0])

        for neighbor in alpha_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 7:
                # Check if the nitrogen is part of -NHOH or -N(OH)2
                
                # Get atoms connected to the nitrogen
                n_neighbors = neighbor.GetNeighbors()

                # Check against pattern 1 or pattern 2
                
                if mol.HasSubstructMatch(n_hydroxy_pattern_1, rootedAtAtom=neighbor.GetIdx()):
                    found_n_hydroxy = True
                    break
                elif mol.HasSubstructMatch(n_hydroxy_pattern_2, rootedAtAtom=neighbor.GetIdx()):
                    found_n_hydroxy = True
                    break
        if found_n_hydroxy:
            break
        
    if not found_n_hydroxy:
        return False, "No N-hydroxy group attached to the alpha amino group found"

    return True, "Contains an alpha amino acid structure with at least one hydroxy group on the amino nitrogen"