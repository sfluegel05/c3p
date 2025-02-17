"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: Lipopeptide – a compound consisting of a peptide with an attached lipid.
A molecule is considered a lipopeptide if it contains at least two amide bonds and at least one 
long aliphatic chain (eight or more contiguous sp3, nonaromatic carbons connected by SINGLE bonds)
directly attached to one of the amide substructures.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def is_valid_chain_atom(atom):
    """Check if an atom is a carbon, nonaromatic and sp3 hybridized."""
    return (atom.GetAtomicNum() == 6 and 
            (not atom.GetIsAromatic()) and 
            atom.GetHybridization().name == "SP3")

def get_chain_length_from_atom(mol, start_idx, visited=None):
    """
    Recursively traverses (using DFS) only along sp3 nonaromatic carbons that are connected by single bonds.
    Returns the maximum number of connected carbons (including the start).
    """
    if visited is None:
        visited = set()
    visited.add(start_idx)
    start_atom = mol.GetAtomWithIdx(start_idx)
    max_length = 1  # count the start atom
    # iterate over neighbors that are linked by a single bond and satisfy valid chain criteria
    for bond in start_atom.GetBonds():
        # only traverse single bonds
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        neighbor = bond.GetOtherAtom(start_atom)
        nidx = neighbor.GetIdx()
        if nidx in visited:
            continue
        if is_valid_chain_atom(neighbor):
            length = 1 + get_chain_length_from_atom(mol, nidx, visited.copy())
            if length > max_length:
                max_length = length
    return max_length

def chain_attached_to_peptide(mol, cutoff=8):
    """
    Checks if at least one long aliphatic chain of cutoff or more sp3 carbon atoms is attached directly
    to any amide bond within the molecule.
    For each matched amide bond, we check the neighbors of the amide nitrogen and the carbonyl carbon.
    From each neighbor that is a valid chain atom (and connected by a SINGLE bond),
    we compute the length of the chain (only walking along sp³ carbons).
    """
    # Define a typical amide fragment (only looking at the N and carbonyl C)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    # For each amide match (tuple of atom indices, usually [N, C, O]), use only N and C
    for match in amide_matches:
        if len(match) < 2:
            continue
        for i in [0,1]:  # check neighbors of the N and the C atom
            atom = mol.GetAtomWithIdx(match[i])
            for bond in atom.GetBonds():
                # Only consider single bonds
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                neighbor = bond.GetOtherAtom(atom)
                # We do not want to walk back into the amide fragment itself.
                if neighbor.GetIdx() in match:
                    continue
                if is_valid_chain_atom(neighbor):
                    chain_length = get_chain_length_from_atom(mol, neighbor.GetIdx())
                    if chain_length >= cutoff:
                        return True
    return False

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    
    A lipopeptide must have a peptide component (at least two amide bonds) and a lipid component –
    a long contiguous aliphatic chain (of at least eight sp³ nonaromatic carbons connected by SINGLE bonds)
    attached directly to one of the amide substructures.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a lipopeptide, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count amide bonds using the standard SMARTS [NX3][CX3](=O)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 2:
        return False, f"Only {len(amide_matches)} amide bond(s) found; need at least two for a peptide component"
    
    # Instead of a global chain search, we now check that at least one amide has a directly attached chain
    if not chain_attached_to_peptide(mol, cutoff=8):
        return False, "No long aliphatic chain (>=8 sp3 carbons) is directly attached to any amide substructure"
    
    # For explanation, count the global (longest) chain as a fallback if desired:
    # (Note: the global longest chain may not be directly attached, so we do not use it as a criterion)
    # We now explain the decision based on the number of amide bonds and presence of an attached chain.
    return True, (f"Contains {len(amide_matches)} amide bonds (peptide component) and "
                  f"a long aliphatic chain (>=8 sp3 carbons) is attached to an amide substructure (lipid component)")

# Example usage (uncomment to run):
# test_smiles = "[H][C@@]1(CCCCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)O1"
# res, reason = is_lipopeptide(test_smiles)
# print(res, reason)