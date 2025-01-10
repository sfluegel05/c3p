"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
"""
from rdkit import Chem
from rdkit.Chem import rdchem

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

    # Define peptide bond pattern (amide bond between C=O and N)
    peptide_bond_pattern = Chem.MolFromSmarts("[$(NC(=O))]")
    if peptide_bond_pattern is None:
        return False, "Failed to create peptide bond SMARTS pattern"

    # Find all peptide bonds
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)

    # Build a graph of peptide bonds to identify connected amino acid residues
    residues = set()
    for match in peptide_bonds:
        n_idx = match[0]
        c_idx = match[1]
        # Add the atoms connected to the nitrogen and carbon to the residues set
        n_atom = mol.GetAtomWithIdx(n_idx)
        c_atom = mol.GetAtomWithIdx(c_idx)
        residues.add(n_atom.GetIdx())
        residues.add(c_atom.GetIdx())

    # Count the number of residues by identifying alpha carbons (C-alpha)
    alpha_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() >= 3:
            has_nitrogen = False
            has_carbonyl = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 7:
                    has_nitrogen = True
                elif neighbor.GetAtomicNum() == 6:
                    for bond in neighbor.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            if bond.GetOtherAtom(neighbor).GetAtomicNum() == 8:
                                has_carbonyl = True
                if has_nitrogen and has_carbonyl:
                    alpha_carbons.append(atom)
                    break

    n_residues = len(alpha_carbons)

    if n_residues != 2:
        return False, f"Found {n_residues} amino acid residues, expected 2 for a dipeptide"

    return True, "Molecule is a dipeptide composed of two amino acid residues connected by peptide linkage"