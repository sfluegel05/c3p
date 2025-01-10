"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: CHEBI:25212 lipopeptide

"""
from rdkit import Chem

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is a compound consisting of a peptide with attached lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify peptide bonds (N-C(=O)-C)
    peptide_bond_pattern = Chem.MolFromSmarts("N[C](=O)[C,N]")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bonds) < 2:
        return False, "Not enough peptide bonds found (need at least 2)"

    # Identify long aliphatic chains (lipid moiety)
    # Define lipid as a chain of at least 8 consecutive non-ring carbons
    aliphatic_chain_pattern = Chem.MolFromSmarts("[C;!R][C;!R]{7,}")
    chains = mol.GetSubstructMatches(aliphatic_chain_pattern)
    if len(chains) == 0:
        return False, "No long aliphatic chain (lipid moiety) found (need at least 8 carbons)"

    # Get indices of peptide atoms
    peptide_atoms = set()
    for match in peptide_bonds:
        peptide_atoms.update(match)

    # Check if any lipid chain is connected to the peptide
    connected = False
    for chain in chains:
        chain_atoms = set(chain)
        for peptide_atom in peptide_atoms:
            for lipid_atom in chain_atoms:
                # Check if there is a path between peptide atom and lipid atom
                try:
                    path = Chem.rdmolops.GetShortestPath(mol, peptide_atom, lipid_atom)
                    if path:
                        connected = True
                        break
                except RuntimeError:
                    continue
            if connected:
                break
        if connected:
            break

    if not connected:
        return False, "Lipid moiety not connected to peptide"

    return True, "Contains peptide bonds and attached lipid moiety"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25212',
        'name': 'lipopeptide',
        'definition': 'A compound consisting of a peptide with attached lipid.',
        'parents': []
    }
}