"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: CHEBI:51728 aliphatic nitrile
An aliphatic nitrile is any nitrile derived from an aliphatic compound. 
Examples include acetonitrile, butyronitrile, and isovaleronitrile.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdchem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for nitrile group
    nitrile_pattern = Chem.MolFromSmarts("[C-]#[N]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if not nitrile_matches:
        return False, "No nitrile group found"
    
    # Check for aliphatic structure
    aromatic_atoms = [atom.GetIsAromatic() for atom in mol.GetAtoms()]
    if any(aromatic_atoms):
        # Molecule contains aromatic atoms, check if nitrile is on aliphatic part
        aliphatic_nitrile = False
        for match in nitrile_matches:
            nitrile_atom = mol.GetAtomWithIdx(match[0])
            is_aliphatic = True
            for neighbor in nitrile_atom.GetNeighbors():
                if neighbor.GetIsAromatic():
                    is_aliphatic = False
                    break
            if is_aliphatic:
                aliphatic_nitrile = True
                break
        if not aliphatic_nitrile:
            return False, "Nitrile group is attached to an aromatic ring"
    
    # Check for common aliphatic nitrile substructures
    aliphatic_patterns = [Chem.MolFromSmarts("[C-](C)(C)C#N"), # isovaleronitrile
                          Chem.MolFromSmarts("CCC#N"), # propionitrile
                          Chem.MolFromSmarts("CC#N"), # acetonitrile
                          Chem.MolFromSmarts("CCCC#N"), # butyronitrile
                          Chem.MolFromSmarts("N[C@@H](CCC#N)C(O)=O")] # alpha-amino-gamma-cyanobutanoic acid
    for pattern in aliphatic_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches common aliphatic nitrile substructure"
    
    # No common substructure found, but molecule contains aliphatic nitrile group
    return True, "Molecule contains an aliphatic nitrile group"