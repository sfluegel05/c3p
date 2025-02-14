"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:38170 quinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is defined as a compound having a fully conjugated cyclic dione structure,
    such as that of benzoquinones, derived from aromatic compounds by conversion of an 
    even number of -CH= groups into -C(=O)- groups with any necessary rearrangement of 
    double bonds (polycyclic and heterocyclic analogues are included).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define quinone pattern
    quinone_pattern = Chem.MolFromSmarts("[#6,#7,#8,#16]1([#6]=[#6])([#8])=[#6][#6]=[#6][#6]1=[#8]")
    
    # Check for quinone pattern
    matches = mol.GetSubstructMatches(quinone_pattern)
    if not matches:
        return False, "No quinone substructure found"
    
    # Check for conjugation and ring membership
    for match in matches:
        is_conjugated = check_conjugation(mol, match)
        is_cyclic = all(mol.GetAtomWithIdx(idx).IsInRing() for idx in match)
        
        if is_conjugated and is_cyclic:
            return True, "Molecule contains a fully conjugated cyclic dione structure (quinone)"
    
    return False, "Quinone substructure is not fully conjugated or not cyclic"

def check_conjugation(mol, atom_indices):
    """
    Checks if a given substructure is fully conjugated in a molecule.
    
    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object
        atom_indices (list): List of atom indices representing the substructure
        
    Returns:
        bool: True if the substructure is fully conjugated, False otherwise
    """
    conjugated_atoms = set(atom_indices)
    for idx in atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in conjugated_atoms and mol.GetBondBetweenAtoms(idx, neighbor.GetIdx()).GetIsConjugated():
                conjugated_atoms.add(neighbor.GetIdx())
    
    return len(conjugated_atoms) == len(atom_indices)