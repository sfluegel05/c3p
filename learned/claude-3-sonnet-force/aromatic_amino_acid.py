"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: CHEBI:33567 aromatic amino acid
An amino acid whose structure includes an aromatic ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for amino acid pattern
    amino_acid_pattern = Chem.MolFromSmarts("[NX3H2][CX4H]([CX3](=O)[OX1-,OX2H1])[CX4H]")
    amino_acid_match = mol.GetSubstructMatches(amino_acid_pattern)
    if not amino_acid_match:
        return False, "No amino acid pattern found"
    
    # Get smallest set of smallest rings (SSSR)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Check for aromatic rings
    aromatic_rings = [ring for ring in ring_info.BondRings() if ring_info.IsBondAromaticRing(ring)]
    if not aromatic_rings:
        return False, "No aromatic rings found"
    
    # Check if aromatic ring is part of the amino acid backbone
    amino_acid_atoms = set.union(*[set(mol.GetAtomWithIdx(idx).GetNeighbors()) for idx in amino_acid_match[0]])
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        if ring_atoms.intersection(amino_acid_atoms):
            return True, "Contains an aromatic ring and an amino acid backbone"
    
    return False, "Aromatic ring not part of amino acid backbone"