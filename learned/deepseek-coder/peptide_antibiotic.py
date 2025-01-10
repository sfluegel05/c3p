"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: CHEBI:36080 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    A peptide antibiotic is a peptide that exhibits antimicrobial properties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bonds (amide bonds, -C(=O)-N-)
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 2:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need at least 2"

    # Count the number of amino acids (approximated by the number of peptide bonds + 1)
    n_amino_acids = len(peptide_bond_matches) + 1
    if n_amino_acids < 4:
        return False, f"Only {n_amino_acids} amino acids, need at least 4"

    # Check for cyclic structure (optional, as some peptide antibiotics are linear)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        # If cyclic, check for at least one large ring (common in cyclic peptide antibiotics)
        large_ring_found = any(len(ring) >= 6 for ring in ring_info.AtomRings())
        if not large_ring_found:
            return False, "No large ring found, cyclic peptide antibiotics typically have large rings"

    # Check for non-standard amino acids (e.g., D-amino acids, modified amino acids)
    # This is a heuristic and may not catch all cases
    non_standard_aa_pattern = Chem.MolFromSmarts("[NX3H][CX4H]([CX4H])[CX3](=[OX1])")
    non_standard_aa_matches = mol.GetSubstructMatches(non_standard_aa_pattern)
    if len(non_standard_aa_matches) == 0:
        return False, "No non-standard amino acids found, peptide antibiotics often contain modified amino acids"

    # Check molecular weight - peptide antibiotics typically have MW > 500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for peptide antibiotic"

    # Count nitrogen and oxygen atoms - peptide antibiotics typically have a high N/O count
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 4:
        return False, "Too few nitrogen atoms for peptide antibiotic"
    if o_count < 4:
        return False, "Too few oxygen atoms for peptide antibiotic"

    return True, "Contains peptide backbone with multiple amino acids, likely a peptide antibiotic"