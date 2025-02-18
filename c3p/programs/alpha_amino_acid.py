"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: CHEBI:33709 alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group and a carboxyl group, with the amino
    group located on the carbon atom adjacent (alpha) to the carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for amino and carboxyl groups
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1&!$(NC=O)]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if not amino_matches or not carboxyl_matches:
        return False, "Missing amino and/or carboxyl group"
    
    # Check if amino group is alpha to carboxyl
    for amino_idx in amino_matches:
        amino_atom = mol.GetAtomWithIdx(amino_idx)
        for bond in amino_atom.GetBonds():
            carbon_atom = bond.GetOtherAtom(amino_atom)
            if carbon_atom.GetDegree() == 3:
                for bond2 in carbon_atom.GetBonds():
                    neighbor = bond2.GetOtherAtom(carbon_atom)
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                        return True, "Contains amino group alpha to carboxyl group"
    
    return False, "Amino group not alpha to carboxyl group"