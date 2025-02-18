"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: CHEBI:17578 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    Must contain:
    1. Sphingosine backbone with amino alcohol and double bond in the chain
    2. Amide-linked fatty acyl group
    3. Minimum 12 carbons in sphingosine chain
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for amide group (N connected to carbonyl)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"

    # Iterate through all amide matches to find valid structure
    for amide_match in amide_matches:
        nitrogen_idx = amide_match[0]
        carbonyl_idx = amide_match[1]
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Find sphingosine carbon (connected to nitrogen and has OH)
        sphingosine_carbon = None
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetIdx() == carbonyl_idx:
                continue  # Skip the carbonyl carbon
            # Check if this neighbor has an OH group
            for nbr_bond in neighbor.GetBonds():
                other_atom = nbr_bond.GetOtherAtom(neighbor)
                if other_atom.GetAtomicNum() == 8 and nbr_bond.GetBondType() == Chem.BondType.SINGLE:
                    sphingosine_carbon = neighbor
                    break
            if sphingosine_carbon:
                break
        if not sphingosine_carbon:
            continue  # No amino alcohol found
        
        # Traverse sphingosine chain to check length and double bond
        visited = set()
        stack = [(sphingosine_carbon, 0)]
        chain_atoms = set()
        double_bond_found = False
        
        while stack:
            atom, depth = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() != 6:
                continue  # Only track carbons
            chain_atoms.add(atom.GetIdx())
            
            # Check bonds for double bonds in the chain
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(atom)
                    if other.GetAtomicNum() == 6 and other.GetIdx() in chain_atoms:
                        double_bond_found = True
            
            # Continue traversal avoiding backtracking to nitrogen/OH
            for bond in atom.GetBonds():
                other = bond.GetOtherAtom(atom)
                if other.GetIdx() == nitrogen_idx:
                    continue  # Don't go back to amide nitrogen
                if other.GetAtomicNum() == 8 and other in [a for a in sphingosine_carbon.GetNeighbors()]:
                    continue  # Skip OH oxygen
                if other.GetAtomicNum() == 6 and other.GetIdx() not in visited:
                    stack.append((other, depth + 1))
        
        # Verify chain requirements
        chain_length = len(chain_atoms)
        if chain_length >= 12 and double_bond_found:
            return True, "Contains sphingosine backbone with N-linked fatty acyl chain"
    
    return False, "Missing sphingosine backbone features or chain too short"