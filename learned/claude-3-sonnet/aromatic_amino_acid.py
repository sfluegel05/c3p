"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid must have an aromatic ring and an amino acid group.

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

    # Check for aromatic ring
    ring_info = mol.GetRingInfo()
    aromatic_rings = []
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings.append(ring)
            
    if not aromatic_rings:
        return False, "No aromatic ring found"

    # Look for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Look for amine group (primary or secondary)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No suitable amine group found"

    # Check for peptide bonds - if present, likely a peptide not an amino acid
    peptide_pattern = Chem.MolFromSmarts("[NX3;H1,H0][CX3](=[OX1])[#6]")
    if len(mol.GetSubstructMatches(peptide_pattern)) > 0:
        return False, "Contains peptide bonds - likely a peptide not an amino acid"

    # Get matches for both amine and carboxyl groups
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)

    # Check if amine and carboxyl are connected through a reasonable path
    # and are properly connected to an aromatic ring
    for amine_match in amine_matches:
        amine_idx = amine_match[0]
        for carboxyl_match in carboxyl_matches:
            carboxyl_idx = carboxyl_match[0]
            path = Chem.GetShortestPath(mol, amine_idx, carboxyl_idx)
            if 2 <= len(path) <= 4:  # Reasonable path length for amino acid
                # Check connection to aromatic ring
                for ring in aromatic_rings:
                    for ring_atom in ring:
                        to_amine = len(Chem.GetShortestPath(mol, ring_atom, amine_idx))
                        to_carboxyl = len(Chem.GetShortestPath(mol, ring_atom, carboxyl_idx))
                        if min(to_amine, to_carboxyl) <= 4:
                            return True, "Contains aromatic ring properly connected to amino acid group"

    return False, "Aromatic ring not properly connected to amino acid group"