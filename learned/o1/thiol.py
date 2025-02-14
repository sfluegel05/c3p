"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: CHEBI:29232 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH) is attached to a carbon atom
    of any aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thiol SMARTS pattern: sulfur atom bonded to hydrogen and carbon
    # [#16]-[H]: sulfur bonded to hydrogen
    # [#16]-[H]-[C]: sulfur bonded to hydrogen and carbon
    thiol_pattern = Chem.MolFromSmarts("[#16][H]")  # Sulfur atom bonded to hydrogen

    # Search for thiol group
    matches = mol.GetSubstructMatches(thiol_pattern)
    if matches:
        for match in matches:
            s_idx = match[0]  # Index of the sulfur atom
            s_atom = mol.GetAtomWithIdx(s_idx)
            # Check if sulfur is bonded to at least one carbon atom
            bonded_to_carbon = False
            for neighbor in s_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon atom
                    bonded_to_carbon = True
                    break
            if bonded_to_carbon:
                # Exclude peptides by checking for peptide bonds (amide linkages)
                # Peptide bond pattern: [C](=O)[N]
                peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
                if mol.HasSubstructMatch(peptide_bond_pattern):
                    return False, "Contains a thiol group but is likely a peptide or protein"
                else:
                    return True, "Contains a thiol group (-SH) attached to a carbon atom"
        return False, "Thiol group (-SH) not attached to a carbon atom"
    else:
        return False, "No thiol group (-SH) found"