"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: CHEBI:24651 organofluorine compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound is defined as a compound containing at least one carbon-fluorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon and fluorine atoms
    has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    has_fluorine = any(atom.GetAtomicNum() == 9 for atom in mol.GetAtoms())
    
    if not has_carbon and not has_fluorine:
        return False, "Contains neither carbon nor fluorine atoms"
    if not has_carbon:
        return False, "Contains fluorine but no carbon atoms"
    if not has_fluorine:
        return False, "Contains no fluorine atoms"

    # Look for any carbon-fluorine bonds (both aromatic and aliphatic)
    cf_bonds = 0
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Check if bond is between C and F
        if ((atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 9) or
            (atom1.GetAtomicNum() == 9 and atom2.GetAtomicNum() == 6)):
            cf_bonds += 1

    if cf_bonds == 0:
        return False, "No carbon-fluorine bonds found"

    # Count total fluorine atoms
    f_atoms = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9])
    
    return True, f"Contains {f_atoms} fluorine atom{'s' if f_atoms>1 else ''} ({cf_bonds} C-F bond{'s' if cf_bonds>1 else ''})"