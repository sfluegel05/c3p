"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: methyl sulfide
"""

from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is any molecule where a sulfur atom is connected via single bonds
    to two carbon atoms, and at least one of those carbons is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Pattern for sulfur connected to a methyl group and any carbon
    methyl_sulfide_pattern = Chem.MolFromSmarts("[CH3]-[S]-[*]")
    # Pattern to identify sulfur connected to two carbons via single bonds
    sulfide_pattern = Chem.MolFromSmarts("[S;D2;$([#6]-S-[#6])]")

    # Search for methyl sulfide pattern
    if mol.HasSubstructMatch(methyl_sulfide_pattern):
        # Check if sulfur atom is connected to exactly two carbons
        sulfide_matches = mol.GetSubstructMatches(sulfide_pattern)
        for match in sulfide_matches:
            sulfur_idx = match[0]  # Index of the sulfur atom
            sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
            neighbors = sulfur_atom.GetNeighbors()
            carbon_neighbors = [atom for atom in neighbors if atom.GetAtomicNum() == 6]

            # Ensure sulfur is connected to exactly two carbons
            if len(carbon_neighbors) != 2:
                continue

            # Check if one of the carbons is a methyl group (degree 1)
            methyl_group_found = False
            for carbon in carbon_neighbors:
                if carbon.GetDegree() == 1:
                    methyl_group_found = True
                    break
            if not methyl_group_found:
                continue

            # Exclude peptides by checking for peptide bonds (amide linkages)
            peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
            if mol.HasSubstructMatch(peptide_bond_pattern):
                # Check if sulfur is part of an amino acid residue (like methionine)
                # Look for C-C-S-C pattern where sulfur is connected to an alpha carbon
                methionine_pattern = Chem.MolFromSmarts("N[C;!R]-[C;!R]-S-C")
                if mol.HasSubstructMatch(methionine_pattern):
                    return False, "Sulfur atom is part of a methionine residue in a peptide"
                else:
                    # Peptide bond present but sulfur not part of amino acid residue
                    return True, "Contains methyl sulfide group"

            # No peptide bonds involving sulfur atom
            return True, "Contains methyl sulfide group"

    return False, "Does not contain the methyl sulfide group as defined"