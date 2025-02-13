"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is a chlorobenzene with three chloro substituents on at least one benzene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify aromatic rings in the molecule
    aromatic_rings = [ring for ring in Chem.GetSymmSSSR(mol) if all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring)]
    
    # Check each aromatic ring for Cl substituents
    for ring in aromatic_rings:
        chloro_count = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check if the atom is a carbon and has a chlorine neighbor
            if atom.GetSymbol() == 'C':
                is_cl_substituent = any(neighbor.GetSymbol() == 'Cl' for neighbor in atom.GetNeighbors())
                if is_cl_substituent:
                    chloro_count += 1
        
        if chloro_count == 3:
            return True, "Contains an aromatic benzene ring with three chloro substituents"
    
    return False, "Benzene ring with three chloro substituents not found"