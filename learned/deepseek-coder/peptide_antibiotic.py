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

    # Check for modified amino acids or unusual structures
    # Look for D-amino acids, N-methylation, or other modifications
    modified_aa_pattern = Chem.MolFromSmarts("[NX3H][CX4H]([CX4H])[CX3](=[OX1])")
    modified_aa_matches = mol.GetSubstructMatches(modified_aa_pattern)
    if len(modified_aa_matches) == 0:
        # Also check for other common modifications
        other_mod_pattern = Chem.MolFromSmarts("[NX3H][CX4H]([CX4H])[CX3](=[OX1])[NX3H]")
        other_mod_matches = mol.GetSubstructMatches(other_mod_pattern)
        if len(other_mod_matches) == 0:
            # Check for D-amino acids
            d_aa_pattern = Chem.MolFromSmarts("[NX3H][CX4H]([CX4H])[CX3](=[OX1])[CX4H]([CX4H])[CX3](=[OX1])")
            d_aa_matches = mol.GetSubstructMatches(d_aa_pattern)
            if len(d_aa_matches) == 0:
                return False, "No modified amino acids found, peptide antibiotics often contain modified amino acids"

    # Check molecular weight - peptide antibiotics typically have MW > 1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for peptide antibiotic"

    # Count nitrogen and oxygen atoms - peptide antibiotics typically have a high N/O count
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 8:
        return False, f"Too few nitrogen atoms ({n_count}) for peptide antibiotic"
    if o_count < 8:
        return False, f"Too few oxygen atoms ({o_count}) for peptide antibiotic"

    # Check for hydrophobic regions (common in peptide antibiotics)
    hydrophobic_pattern = Chem.MolFromSmarts("[CX4H][CX4H][CX4H][CX4H]")
    hydrophobic_matches = mol.GetSubstructMatches(hydrophobic_pattern)
    if len(hydrophobic_matches) < 2:
        return False, "Insufficient hydrophobic regions for peptide antibiotic"

    # Check for cationic residues (lysine, arginine)
    cationic_pattern = Chem.MolFromSmarts("[NX3H2][CX4H2][CX4H2][CX4H2][NX3H2]")
    cationic_matches = mol.GetSubstructMatches(cationic_pattern)
    if len(cationic_matches) == 0:
        return False, "No cationic residues found, peptide antibiotics often contain lysine or arginine"

    # Check for heterocycles (thiazole, oxazole, etc.)
    heterocycle_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#7][#6]1")
    heterocycle_matches = mol.GetSubstructMatches(heterocycle_pattern)
    if len(heterocycle_matches) == 0:
        return False, "No heterocycles found, peptide antibiotics often contain thiazole or oxazole rings"

    return True, "Contains peptide backbone with multiple amino acids, modifications, and structural features typical of peptide antibiotics"