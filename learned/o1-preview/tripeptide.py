"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is an oligopeptide that consists of three amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tripeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Function to identify backbone peptide bonds
    def get_peptide_bonds(mol):
        peptide_bond_idxs = []
        for bond in mol.GetBonds():
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            # Look for C-N single bonds
            if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 7) or \
               (begin_atom.GetAtomicNum() == 7 and end_atom.GetAtomicNum() == 6):
                C_atom = begin_atom if begin_atom.GetAtomicNum() == 6 else end_atom
                N_atom = begin_atom if begin_atom.GetAtomicNum() == 7 else end_atom
                # Check if C atom is carbonyl carbon (connected to O via double bond)
                is_carbonyl = False
                for neighbor in C_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8:
                        bond_to_O = mol.GetBondBetweenAtoms(C_atom.GetIdx(), neighbor.GetIdx())
                        if bond_to_O.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            is_carbonyl = True
                            break
                if not is_carbonyl:
                    continue
                # Check that N atom is connected to alpha carbon (C)
                N_connected_to_Calpha = False
                for neighbor in N_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != C_atom.GetIdx():
                        N_connected_to_Calpha = True
                        break
                if not N_connected_to_Calpha:
                    continue
                # Check that C atom is connected to alpha carbon (C)
                C_connected_to_Calpha = False
                for neighbor in C_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != N_atom.GetIdx():
                        if neighbor.GetAtomicNum() == 6:
                            C_connected_to_Calpha = True
                            break
                if not C_connected_to_Calpha:
                    continue
                peptide_bond_idxs.append(bond.GetIdx())
        return peptide_bond_idxs

    peptide_bonds = get_peptide_bonds(mol)

    if len(peptide_bonds) < 2:
        return False, f"Found {len(peptide_bonds)} peptide bonds, expected at least 2 for tripeptide"

    # Break molecule at peptide bonds
    fragments = Chem.FragmentOnBonds(mol, peptide_bonds, addDummies=True)
    frags = Chem.GetMolFrags(fragments, asMols=True)
    num_fragments = len(frags)

    # The number of amino acid residues is number of peptide bonds + 1
    num_residues = len(peptide_bonds) + 1

    if num_residues != 3:
        return False, f"Found {num_residues} amino acid residues, expected 3 for tripeptide"

    return True, "Molecule is a tripeptide with three amino acid residues connected by peptide bonds"