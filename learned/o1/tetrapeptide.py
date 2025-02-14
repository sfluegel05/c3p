"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is a molecule containing four amino acid residues connected by peptide linkages.

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

    try:
        # Define peptide bond pattern (amide bond between N and C=O)
        peptide_bond_pattern = Chem.MolFromSmarts("[$(NC(=O))]")

        # Find all peptide bonds in the molecule
        peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)

        # Count the number of peptide bonds
        num_peptide_bonds = len(peptide_bonds)

        if num_peptide_bonds != 3:
            return False, f"Found {num_peptide_bonds} peptide bonds, expected 3 for tetrapeptide"

        # Optional: Ensure peptide bonds form a continuous chain
        # Check that peptide bonds are connected in sequence
        amide_bond_idxs = []
        for match in peptide_bonds:
            n_idx = match[0]
            c_idx = match[1]

            # Find the bond between N and C=O (amide bond)
            bond = mol.GetBondBetweenAtoms(n_idx, c_idx)
            if bond is not None:
                amide_bond_idxs.append(bond.GetIdx())

        # Get the subgraph formed by peptide bonds
        fragmented_mol = Chem.FragmentOnBonds(mol, amide_bond_idxs, addDummies=False)
        frags = Chem.GetMolFrags(fragmented_mol, asMols=True)

        # Count the number of fragments (should be equal to number of residues)
        num_fragments = len(frags)

        if num_fragments != 4:
            return False, f"Found {num_fragments} amino acid residues, expected 4 for tetrapeptide"

        return True, "Molecule is a tetrapeptide with 4 amino acid residues connected via peptide bonds"

    except Exception as e:
        return False, f"Error during analysis: {str(e)}"