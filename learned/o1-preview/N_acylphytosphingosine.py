"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: CHEBI:76955 N-acylphytosphingosine
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide that is phytosphingosine
    having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
        
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for amide bond pattern (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found"
    
    # Loop through all amide bonds to check for correct connections
    for match in amide_matches:
        carbonyl_c_idx = match[0]
        nitrogen_idx = match[2]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)

        # Get neighbors of nitrogen atom excluding the carbonyl carbon
        neighbors = [a for a in nitrogen_atom.GetNeighbors() if a.GetIdx() != carbonyl_c_idx]
        if len(neighbors) == 0:
            continue  # No other atoms attached to nitrogen

        # Assume the first neighbor is alpha carbon
        alpha_c_atom = neighbors[0]

        # Check for hydroxyl on alpha carbon
        has_alpha_oh = False
        for neighbor in alpha_c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                bond = mol.GetBondBetweenAtoms(alpha_c_atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    has_alpha_oh = True
                    break
        if not has_alpha_oh:
            continue  # No hydroxyl on alpha carbon

        # Get beta carbon(s)
        beta_c_atoms = [a for a in alpha_c_atom.GetNeighbors() if a.GetIdx() != nitrogen_idx and a.GetAtomicNum() == 6]
        if len(beta_c_atoms) == 0:
            continue  # No beta carbon
        beta_c_atom = beta_c_atoms[0]

        # Check for hydroxyl on beta carbon
        has_beta_oh = False
        for neighbor in beta_c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() != alpha_c_atom.GetIdx():
                bond = mol.GetBondBetweenAtoms(beta_c_atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    has_beta_oh = True
                    break
        if not has_beta_oh:
            continue  # No hydroxyl on beta carbon

        # Get gamma carbon(s)
        gamma_c_atoms = [a for a in beta_c_atom.GetNeighbors() if a.GetIdx() != alpha_c_atom.GetIdx() and a.GetAtomicNum() == 6]
        if len(gamma_c_atoms) == 0:
            continue  # No gamma carbon
        gamma_c_atom = gamma_c_atoms[0]

        # Check for hydroxyl on gamma carbon
        has_gamma_oh = False
        for neighbor in gamma_c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() != beta_c_atom.GetIdx():
                bond = mol.GetBondBetweenAtoms(gamma_c_atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    has_gamma_oh = True
                    break
        if not has_gamma_oh:
            continue  # No hydroxyl on gamma carbon

        # From gamma carbon, check for long carbon chain (minimum 10 carbons)
        chain_length = 0
        visited = set([nitrogen_idx, alpha_c_atom.GetIdx(), beta_c_atom.GetIdx(), gamma_c_atom.GetIdx()])
        current_atom = gamma_c_atom
        while True:
            neighbors = [a for a in current_atom.GetNeighbors() if a.GetIdx() not in visited and a.GetAtomicNum() == 6]
            if len(neighbors) == 0:
                break
            chain_length += 1
            current_atom = neighbors[0]
            visited.add(current_atom.GetIdx())

        if chain_length < 10:
            continue  # Aliphatic chain too short

        # If all conditions are met
        return True, "Molecule is an N-acylphytosphingosine"

    return False, "No matching phytosphingosine backbone attached via amide bond"