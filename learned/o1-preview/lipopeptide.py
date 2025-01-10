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

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify peptide bonds (amide bonds connecting amino acids)
    peptide_bond_pattern = Chem.MolFromSmarts("N[C](=O)C")  # N-C(=O)-C
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bonds) < 2:
        return False, "Not enough peptide bonds found (need at least 2)"
    
    # Identify aliphatic carbons (including in rings and chains)
    aliphatic_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not aliphatic_carbons:
        return False, "No aliphatic carbons found"
    
    # Check if any aliphatic carbon is connected to the peptide
    peptide_atoms = set()
    for match in peptide_bonds:
        peptide_atoms.update(match)
    
    connected = False
    for carbon_atom in aliphatic_carbons:
        for peptide_atom in peptide_atoms:
            if carbon_atom != peptide_atom:
                path = Chem.rdmolops.GetShortestPath(mol, carbon_atom, peptide_atom)
                if path and len(path) > 0:
                    connected = True
                    break
        if connected:
            break
    
    if not connected:
        return False, "No lipid moiety connected to peptide found"
    
    return True, "Contains peptide bonds and attached lipid moiety"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25212',
        'name': 'lipopeptide',
        'definition': 'A compound consisting of a peptide with attached lipid.',
        'parents': []
    }
}