"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: N-acetyl amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl amino acid based on its SMILES string.
    An N-acetyl amino acid is an amino acid where an acetyl group is attached to the nitrogen atom of the amino group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acetyl amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # SMARTS pattern for N-acetyl group attached to nitrogen
    # Matches nitrogen attached to a carbonyl carbon (C=O) bonded to a methyl group
    n_acetyl_pattern = Chem.MolFromSmarts("N[C](=O)C")
    acetyl_matches = mol.GetSubstructMatches(n_acetyl_pattern)
    if not acetyl_matches:
        return False, "No N-acetyl group found"

    # SMARTS pattern for amino acid backbone
    # Matches alpha carbon connected to a carboxyl group (C(=O)O)
    amino_acid_backbone_pattern = Chem.MolFromSmarts("[CX4H1,CX4H2][CX3](=O)[O;H1,H0-]")
    backbone_matches = mol.GetSubstructMatches(amino_acid_backbone_pattern)
    if not backbone_matches:
        return False, "No amino acid backbone found"

    # Check that the acetylated nitrogen is part of the amino acid backbone
    for n_acetyl_match in acetyl_matches:
        acetyl_nitrogen_idx = n_acetyl_match[0]
        # Find alpha carbon attached to the nitrogen
        neighbor_atoms = mol.GetAtomWithIdx(acetyl_nitrogen_idx).GetNeighbors()
        for atom in neighbor_atoms:
            if atom.GetAtomicNum() == 6:  # Carbon atom
                alpha_carbon_idx = atom.GetIdx()
                # Check if this carbon is part of the amino acid backbone
                for backbone_match in backbone_matches:
                    if alpha_carbon_idx == backbone_match[0]:
                        return True, "Molecule is an N-acetyl amino acid"
    return False, "Acetyl group not attached to amino nitrogen of amino acid backbone"