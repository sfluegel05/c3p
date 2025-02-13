"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: isothiocyanate compounds (R-N=C=S)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    Isothiocyanates have the general formula R-N=C=S where R is an organic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N=C=S group in both possible resonance forms
    isothiocyanate_pattern1 = Chem.MolFromSmarts("[NX2]=[CX2]=[SX1]")  # N=C=S form
    isothiocyanate_pattern2 = Chem.MolFromSmarts("[NX2-][CX2]#[SX1+]") # N-C≡S+ form
    
    matches = []
    for pattern in [isothiocyanate_pattern1, isothiocyanate_pattern2]:
        if pattern is not None:
            matches.extend(mol.GetSubstructMatches(pattern))
    
    if not matches:
        return False, "No isothiocyanate group found"
    
    valid_groups = 0
    for match in matches:
        n_atom, c_atom, s_atom = [mol.GetAtomWithIdx(idx) for idx in match]
        
        # Detailed validation of the isothiocyanate group
        valid = True
        
        # Check nitrogen
        if n_atom.GetDegree() != 2 or n_atom.GetTotalValence() not in [3, 4]:
            valid = False
            continue
            
        # Check carbon
        if c_atom.GetDegree() != 2 or c_atom.GetTotalValence() not in [4]:
            valid = False
            continue
            
        # Check sulfur
        if s_atom.GetDegree() != 1 or s_atom.GetTotalValence() not in [2]:
            valid = False
            continue
            
        # Verify R group is organic (connected to nitrogen)
        r_atoms = [atom for atom in n_atom.GetNeighbors() if atom.GetIdx() != c_atom.GetIdx()]
        if not r_atoms:
            valid = False
            continue
            
        r_atom = r_atoms[0]
        
        # R group should be carbon or connected to carbon
        is_organic = False
        if r_atom.GetAtomicNum() == 6:  # Carbon
            is_organic = True
        else:
            # Check if R is connected to any carbon
            for neighbor in r_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    is_organic = True
                    break
        
        if not is_organic:
            valid = False
            continue
            
        # Check bond types
        n_c_bond = mol.GetBondBetweenAtoms(n_atom.GetIdx(), c_atom.GetIdx())
        c_s_bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), s_atom.GetIdx())
        
        if n_c_bond is None or c_s_bond is None:
            valid = False
            continue
            
        # Valid bond types combinations:
        # N=C=S or N-C≡S
        valid_bonds = (
            (n_c_bond.GetBondType() == Chem.BondType.DOUBLE and c_s_bond.GetBondType() == Chem.BondType.DOUBLE) or
            (n_c_bond.GetBondType() == Chem.BondType.SINGLE and c_s_bond.GetBondType() == Chem.BondType.TRIPLE)
        )
        
        if not valid_bonds:
            valid = False
            continue
            
        if valid:
            valid_groups += 1
    
    if valid_groups == 0:
        return False, "No valid isothiocyanate groups found"
        
    if valid_groups == 1:
        return True, "Contains one isothiocyanate group (R-N=C=S)"
    else:
        return True, f"Contains {valid_groups} isothiocyanate groups"