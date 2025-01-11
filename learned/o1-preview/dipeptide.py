"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is any molecule that contains two amino-acid residues connected by peptide linkages.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify peptide bonds (amide bonds connecting amino acids)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    matches = mol.GetSubstructMatches(peptide_bond_pattern)
    n_peptide_bonds = len(matches)
    if n_peptide_bonds == 0:
        return False, "No peptide bonds found"
    elif n_peptide_bonds > 1:
        return False, f"Found {n_peptide_bonds} peptide bonds, expected 1 for a dipeptide"

    # Fragment the molecule at peptide bonds to get residues
    bonds_to_break = []
    for match in matches:
        c_idx = match[0]
        n_idx = match[1]
        bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
        if bond is not None:
            bonds_to_break.append(bond.GetIdx())
    # Fragment the molecule at the peptide bond
    mol_frag = Chem.FragmentOnBonds(mol, bonds_to_break, addDummies=True)
    # Get the fragments (residues)
    frags = Chem.GetMolFrags(mol_frag, asMols=True)
    n_residues = len(frags)
    if n_residues != 2:
        return False, f"Found {n_residues} residues after fragmentation, expected 2"

    # Check if each fragment is an amino acid residue
    # Simplified check: fragment should have amino and carboxyl groups
    amino_pattern = Chem.MolFromSmarts("[NX3H2,NX3H][CX4]")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    for i, frag in enumerate(frags):
        has_amino = frag.HasSubstructMatch(amino_pattern)
        has_carboxyl = frag.HasSubstructMatch(carboxyl_pattern)
        if not (has_amino and has_carboxyl):
            return False, f"Fragment {i+1} is not an amino acid residue"
    return True, "Molecule is a dipeptide composed of two amino acid residues connected by a peptide bond"