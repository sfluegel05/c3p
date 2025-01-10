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

    # Look for carbon atoms bonded to fluorine atoms
    cf_bond_pattern = Chem.MolFromSmarts("[C]-[F]")
    cf_bond_matches = mol.GetSubstructMatches(cf_bond_pattern)
    
    # Look for CF2 and CF3 groups which may be written differently in SMILES
    cf2_pattern = Chem.MolFromSmarts("[C](F)(F)")
    cf3_pattern = Chem.MolFromSmarts("[C](F)(F)F")
    
    cf2_matches = mol.GetSubstructMatches(cf2_pattern)
    cf3_matches = mol.GetSubstructMatches(cf3_pattern)
    
    total_f_bonds = len(cf_bond_matches) + 2*len(cf2_matches) + 3*len(cf3_matches)
    
    if total_f_bonds == 0:
        # Check if molecule contains both C and F atoms but no C-F bonds
        has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
        has_fluorine = any(atom.GetAtomicNum() == 9 for atom in mol.GetAtoms())
        
        if has_carbon and has_fluorine:
            return False, "Contains C and F atoms but no C-F bonds"
        elif has_carbon and not has_fluorine:
            return False, "Contains no fluorine atoms"
        elif not has_carbon and has_fluorine:
            return False, "Contains fluorine but no carbon atoms"
        else:
            return False, "Contains neither carbon nor fluorine atoms"
    
    # Count number of fluorine atoms and their connections
    f_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
    
    details = []
    if len(cf_bond_matches) > 0:
        details.append(f"{len(cf_bond_matches)} C-F bond{'s' if len(cf_bond_matches)>1 else ''}")
    if len(cf2_matches) > 0:
        details.append(f"{len(cf2_matches)} CF2 group{'s' if len(cf2_matches)>1 else ''}")
    if len(cf3_matches) > 0:
        details.append(f"{len(cf3_matches)} CF3 group{'s' if len(cf3_matches)>1 else ''}")
        
    return True, f"Contains {f_atoms} fluorine atom{'s' if f_atoms>1 else ''} ({', '.join(details)})"