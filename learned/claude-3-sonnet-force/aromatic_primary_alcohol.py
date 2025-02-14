"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: CHEBI:51831 aromatic primary alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol is a primary alcohol in which the alcoholic hydroxy group
    is attached to a carbon which is itself bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove explicit hydrogens without sanitizing (to handle tautomers)
    mol = Chem.RemoveHs(mol, sanitizeRemoving=False)
    
    # Get aromatic atoms
    aromatic_atoms = Chem.GetAromaticAtoms(mol)
    if not aromatic_atoms:
        return False, "No aromatic rings found"
    
    # Get primary alcohol groups
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2][OX2H]")
    primary_alcohol_matches = mol.GetSubstructMatches(primary_alcohol_pattern)
    
    # Check if any primary alcohol group is attached to an aromatic ring
    for match in primary_alcohol_matches:
        alcohol_carbon_idx = match[0]
        alcohol_oxygen_idx = match[1]
        alcohol_carbon = mol.GetAtomWithIdx(alcohol_carbon_idx)
        for neighbor_idx in alcohol_carbon.GetNeighbors():
            neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
            if neighbor_atom.GetIdx() in aromatic_atoms:
                return True, "Primary alcohol group attached to an aromatic ring"
    
    return False, "No primary alcohol group attached to an aromatic ring"