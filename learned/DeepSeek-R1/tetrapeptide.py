"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide contains four amino acid residues connected by three or four peptide bonds (if cyclic).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern to match peptide bonds (amide groups in backbone)
    # This matches [C]=O-N where N is connected to next residue's alpha carbon
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])-[NX3]")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bonds)
    
    # Linear tetrapeptides have 3 peptide bonds, cyclic have 4
    if num_peptide_bonds == 3 or num_peptide_bonds == 4:
        return True, f"Found {num_peptide_bonds} peptide bonds"
    else:
        return False, f"Found {num_peptide_bonds} peptide bonds, need 3 or 4"