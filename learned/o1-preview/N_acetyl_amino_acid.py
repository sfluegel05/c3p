"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is an amino acid where the amino group is acetylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Pattern for acetylated nitrogen connected to alpha carbon
    acetylated_nitrogen_pattern = Chem.MolFromSmarts('N(C(=O)C)-[C]')
    # Pattern for carboxyl group connected to alpha carbon
    carboxyl_group_pattern = Chem.MolFromSmarts('[C](=O)[O-,OH]')

    # Find all matches for acetylated nitrogen connected to alpha carbon
    matches = mol.GetSubstructMatches(acetylated_nitrogen_pattern)

    if not matches:
        return False, "No acetylated amino group found"

    for match in matches:
        nitrogen_idx = match[0]
        alpha_carbon_idx = match[2]
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)

        # Check if alpha carbon is connected to a carboxyl group
        has_carboxyl = False
        for neighbor in alpha_carbon.GetNeighbors():
            if neighbor.GetIdx() != nitrogen_idx:
                if neighbor.GetAtomicNum() == 6:
                    # Check if neighbor is carbonyl carbon
                    for bond in neighbor.GetBonds():
                        other_atom = bond.GetOtherAtom(neighbor)
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and other_atom.GetAtomicNum() == 8:
                            # Found C=O
                            # Now check for attachment to O or [O-]
                            for oxy_neighbor in neighbor.GetNeighbors():
                                if oxy_neighbor.GetAtomicNum() == 8 and oxy_neighbor.GetIdx() != other_atom.GetIdx():
                                    has_carboxyl = True
                                    break
                        if has_carboxyl:
                            break
                elif neighbor.GetAtomicNum() == 8:
                    # Check if it's a carboxylate oxygen directly connected
                    has_carboxyl = True
                    break
            if has_carboxyl:
                break
        
        if has_carboxyl:
            return True, "Molecule matches N-acetyl-amino acid pattern"
    
    return False, "Molecule does not match N-acetyl-amino acid pattern"