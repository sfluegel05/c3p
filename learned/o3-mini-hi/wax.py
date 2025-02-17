"""
Classifies: CHEBI:73702 wax
"""
#!/usr/bin/env python3
"""
Classifies: Wax
Definition: A chemical substance that is an organic compound or mixture of compounds that is composed of long‐chain molecules and is malleable at ambient temperatures.
Approach:
 - Look for at least one ester functional group (–C(=O)O–).
 - For each ester found, identify the two sides:
     • The acyl (acid) part: the atom attached to the carbonyl carbon (other than the ester oxygen).
     • The alkoxy (alcohol) part: the atom attached to the ester oxygen (other than the carbonyl).
 - For each side, perform a depth-first search (considering only carbon atoms) to find the longest contiguous carbon chain.
 - Qualify as a wax if one ester is found where both chains are long (e.g. ≥8 carbons), the overall molecular weight is high (≥300 Da) and the molecule has a good number of rotatable bonds.
 
If any of these criteria fail, the molecule is not classified as a wax.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol, atom_idx, visited=None):
    """
    Recursively calculates the length of the longest contiguous chain of carbon atoms
    starting from the atom with index atom_idx. Only atoms with atomic number 6 (C) are included.
    
    Args:
        mol (Chem.Mol): The RDKit molecule.
        atom_idx (int): Index of the starting atom.
        visited (set): Set of atom indices that have been visited.
    
    Returns:
        int: The length (number of carbon atoms) in the longest chain from the starting atom.
    """
    if visited is None:
        visited = set()
    visited.add(atom_idx)
    atom = mol.GetAtomWithIdx(atom_idx)
    max_length = 1  # count the starting carbon

    for neighbor in atom.GetNeighbors():
        n_idx = neighbor.GetIdx()
        if neighbor.GetAtomicNum() == 6 and n_idx not in visited:
            # Pass a copy of visited to avoid sharing across branches
            branch_length = 1 + longest_carbon_chain(mol, n_idx, visited.copy())
            if branch_length > max_length:
                max_length = branch_length
    return max_length

def is_wax(smiles: str):
    """
    Determines if a molecule qualifies as a wax.
    
    The heuristic criteria used are:
      - The molecule must be organic and parsed correctly.
      - It should contain at least one ester group, defined as a carbonyl directly bound to an oxygen (i.e. [CX3](=O)[OX2]).
      - For at least one ester group, both the acyl (fatty acid) and alkoxy (fatty alcohol) parts must have long carbon chains (here we require at least 8 carbons in their maximum contiguous chain).
      - The overall molecular weight should be on the higher side (here ≥300 Da).
      - The molecule should have a good number of rotatable bonds (≥5) to be consistent with long, flexible chains.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a wax, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Must contain carbon
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not an organic molecule (no carbon atoms)."
    
    # Check molecular weight; waxes tend to have high molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a typical wax."

    # Check the number of rotatable bonds; long chains give a higher count.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, f"Too few rotatable bonds ({n_rotatable}) for a long-chain wax."

    # Define an ester pattern that captures a carbonyl (C=O) attached to an oxygen (excluding acids which bear H)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found in molecule, which is required for a wax."

    # For each ester match, attempt to verify that both sides have long carbon chains.
    for match in ester_matches:
        # match[0]: carbonyl carbon, match[1]: the oxygen of the ester linkage
        carbonyl_idx = match[0]
        oxygen_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)

        # For the acyl (fatty acid) side:
        # Get the neighbor of the carbonyl carbon that is NOT the oxygen in the ester motif.
        acid_chain_length = 0
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetIdx() != oxygen_idx and neighbor.GetAtomicNum() == 6:
                acid_chain_length = longest_carbon_chain(mol, neighbor.GetIdx())
                break  # only consider one acyl chain

        # For the alcohol (fatty alcohol) side:
        alcohol_chain_length = 0
        for neighbor in oxygen_atom.GetNeighbors():
            if neighbor.GetIdx() != carbonyl_idx and neighbor.GetAtomicNum() == 6:
                alcohol_chain_length = longest_carbon_chain(mol, neighbor.GetIdx())
                break

        # Check if both sides have a sufficiently long carbon chain.
        if acid_chain_length >= 8 and alcohol_chain_length >= 8:
            return True, (f"Found an ester group with acyl chain length {acid_chain_length} and "
                          f"alcohol chain length {alcohol_chain_length}, consistent with a wax.")
    return False, "No ester group found with two long-chain carbon moieties."

# Example usage:
if __name__ == "__main__":
    # Example: 2-palmitoyloxypalmityl palmitate (one of the wax examples)
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OCC(CCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
    classification, reason = is_wax(test_smiles)
    print(f"Classification: {classification}\nReason: {reason}")