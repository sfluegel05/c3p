"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: Aromatic amino acids (CHEBI:47769)
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid has an amino group, a carboxylic acid group, and an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for amino group (-NH2) and carboxylic acid (-COOH)
    amino_pattern = Chem.MolFromSmarts("[NH2]")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group"
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group"
    
    # Check amino acid backbone: NH2-C-COOH with at least one additional substituent
    backbone_pattern = Chem.MolFromSmarts("[NH2][C]([CX4])([C]=O)O")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "Not a standard amino acid backbone"
    
    # Check for aromatic rings
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() 
                      if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)]
    if not aromatic_rings:
        return False, "No aromatic rings detected"
    
    return True, "Contains amino acid backbone and aromatic ring(s)"