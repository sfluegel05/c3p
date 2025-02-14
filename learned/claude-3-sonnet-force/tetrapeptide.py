"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: CHEBI:25104 tetrapeptide
A tetrapeptide is any molecule that contains four amino-acid residues connected by peptide linkages.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Define SMARTS patterns for common proteinogenic amino acids
amino_acid_patterns = [
    "[NX3H2,NX4H3+][C@H](C(=O)O)[C@@H](C)O",  # Ala
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(C)=O",  # Val
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(C)(C)C",  # Leu
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(C)(C)O",  # Ile
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)O",  # Asp
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N",  # Asn
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](C(=O)O)CCC(C)C",  # Lys
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](CC(C)C)C(=O)O",  # Arg
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](CC[C@@H](C(=O)O)N)C(=O)O",  # Glu
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](CC(=O)N)C(=O)O",  # Gln
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(C)(C)S",  # Met
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](CC[C@@H](C(=O)O)O)C(=O)O",  # Ser
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](CO)C(=O)O",  # Thr
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](Cc1c[nH]c2c1cccc2)C(=O)O",  # Trp
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](Cc1ccc(cc1)O)C(=O)O",  # Tyr
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](Cc1cnc[nH]1)C(=O)O",  # His
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(C)(C)O",  # Thr
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](CS)C(=O)O",  # Cys
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(C)(C)N",  # Pro
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](CCC(=O)O)C(=O)O",  # Glu
    "[NX3H2,NX4H3+][C@H](C(=O)O)[CH2]C(=O)N[C@@H](CCC(=O)N)C(=O)O"  # Gln
]

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is any molecule that contains four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 3 for tetrapeptide"

    # Check for amino acid residues
    residues = []
    for pattern in amino_acid_patterns:
        residue_pattern = Chem.MolFromSmarts(pattern)
        residue_matches = mol.GetSubstructMatches(residue_pattern)
        residues.extend(residue_matches)

    if len(set(residues)) != 4:
        return False, f"Found {len(set(residues))} amino acid residues, need exactly 4 for tetrapeptide"

    # Check for correct N-terminus and C-terminus
    n_terminus_pattern = Chem.MolFromSmarts("[NX3H2]")
    c_terminus_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    n_terminus_matches = mol.GetSubstructMatches(n_terminus_pattern)
    c_terminus_matches = mol.GetSubstructMatches(c_terminus_pattern)

    if len(n_terminus_matches) != 1 or len(c_terminus_matches) != 1:
        return False, "Incorrect N-terminus or C-terminus"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 800:
        return False, "Molecular weight outside typical range for tetrapeptides"

    return True, "Molecule contains four amino acid residues connected by peptide linkages"