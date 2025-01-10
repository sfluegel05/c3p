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

    # Look for both resonance forms of isothiocyanate group
    patterns = [
        Chem.MolFromSmarts("[NX2]=[CX2]=[SX1]"),  # N=C=S form
        Chem.MolFromSmarts("[NX2]=[CX2]-[SX2-]"),  # Resonance form
    ]
    
    all_matches = []
    for pattern in patterns:
        if pattern is not None:
            matches = mol.GetSubstructMatches(pattern)
            all_matches.extend(matches)
    
    if not all_matches:
        return False, "No isothiocyanate group found"
    
    valid_groups = 0
    for match in all_matches:
        n_atom, c_atom, s_atom = [mol.GetAtomWithIdx(idx) for idx in match]
        
        # Get R group (atom connected to nitrogen that's not the central carbon)
        r_atoms = [atom for atom in n_atom.GetNeighbors() if atom.GetIdx() != c_atom.GetIdx()]
        if not r_atoms:
            continue
        r_atom = r_atoms[0]
        
        # Validate the R group - must be carbon in an organic group
        if r_atom.GetAtomicNum() != 6:  # Must be carbon
            continue
            
        # Check that R is part of an organic group (has other carbons or reasonable substituents)
        r_neighbors = r_atom.GetNeighbors()
        has_carbon_chain = False
        for neighbor in r_neighbors:
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != c_atom.GetIdx():
                has_carbon_chain = True
                break
        if not has_carbon_chain and len(r_neighbors) < 2:
            continue
        
        # Validate the N=C=S connectivity
        n_bonds = [bond.GetBondType() for bond in n_atom.GetBonds()]
        c_bonds = [bond.GetBondType() for bond in c_atom.GetBonds()]
        s_bonds = [bond.GetBondType() for bond in s_atom.GetBonds()]
        
        # Check correct bonding pattern
        if (Chem.BondType.DOUBLE not in n_bonds or 
            len(c_bonds) != 2 or 
            len(s_bonds) > 1):
            continue
            
        valid_groups += 1
    
    if valid_groups == 0:
        return False, "No valid isothiocyanate groups found"
        
    # Return appropriate message based on number of groups
    if valid_groups == 1:
        return True, "Contains one isothiocyanate group (R-N=C=S)"
    else:
        return True, f"Contains {valid_groups} isothiocyanate groups"