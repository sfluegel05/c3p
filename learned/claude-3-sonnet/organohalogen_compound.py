"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: CHEBI:33383 organohalogen compound
An organohalogen compound is defined as a compound containing at least one carbon-halogen bond
(where X is a halogen atom), and the halogen atom is part of the parent molecule.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if molecule contains any halogen atoms as part of the parent molecule
    halogen_atoms = ['Cl', 'Br', 'I', 'F']
    parent_mol_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomMapNum() == 0]
    has_halogen = any(atom.GetSymbol() in halogen_atoms for atom in parent_mol_atoms)
    if not has_halogen:
        return False, "No halogen atoms present in the parent molecule"
    
    # Check if any halogen atom is bonded to a carbon in the parent molecule
    has_c_x_bond = any(atom.GetSymbol() in halogen_atoms and
                       any(nbr_atom.GetSymbol() == 'C' and nbr_atom.GetAtomMapNum() == 0 for nbr_atom in atom.GetNeighbors())
                       for atom in parent_mol_atoms)
    
    if has_c_x_bond:
        return True, "Contains at least one carbon-halogen bond in the parent molecule"
    else:
        return False, "No carbon-halogen bonds found in the parent molecule"