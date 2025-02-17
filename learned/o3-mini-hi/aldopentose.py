"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: aldopentose, defined as 'A pentose with a (potential) aldehyde group at one end.'

An aldopentose can appear as an open‐chain (with a free aldehyde group) or as a cyclic sugar
(a hemiacetal) that is in equilibrium with its open‐chain form. To be an aldopentose the molecule’s 
“main‐chain” should contain exactly 5 carbons and 5 oxygens (as for C5H10O5). 

To avoid false positives (e.g. lactones) we:
  • compute the longest contiguous carbon chain (only carbons connected by C–C bonds) – 
    this should be 5 for a pentose.
  • require the overall count of oxygen atoms to be 5.
  • if the open‐chain aldehyde (SMARTS “[CX3H1](=O)”) is detected, the sugar is classified as aldopentose;
    if not, then (provided the molecule is cyclic and does not contain a lactone (ester C=O–O pattern)),
    it is assumed to be a cyclized aldopentose.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the longest contiguous path (chain) through carbon atoms only.
    In the sugar backbone, the longest carbon chain should be 5 for an aldopentose.
    """
    # Get the indices of all carbon atoms.
    carbon_inds = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Build a graph (dictionary) where each key is a carbon atom index and 
    # each value is a list of neighboring carbon atom indices.
    graph = {idx: [] for idx in carbon_inds}
    for idx in carbon_inds:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                graph[idx].append(nbr.GetIdx())
    max_length = 0
    # A simple DFS to go through all paths (since molecules are small this is acceptable)
    def dfs(current, visited):
        nonlocal max_length
        # Update the maximum path length encountered so far.
        if len(visited) > max_length:
            max_length = len(visited)
        for neighbor in graph[current]:
            if neighbor not in visited:
                dfs(neighbor, visited | {neighbor})
    for idx in carbon_inds:
        dfs(idx, {idx})
    return max_length

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    
    An aldopentose is defined as a 5‐carbon sugar with a (potential) aldehyde group at one end.
    It can exist in an open-chain form, in which an aldehyde is detected, or in a cyclic
    (hemiacetal) form which is in equilibrium with the open-chain form. 

    The classification uses several criteria:
      1. The longest contiguous carbon chain (using only C–C bonds) is expected to have exactly 5 carbons.
      2. The molecule must have exactly 5 oxygen atoms (consistent with C5H10O5).
      3. If an open-chain aldehyde group (SMARTS "[CX3H1](=O)") is found on the structure, 
         the molecule is classified as an open‐chain aldopentose.
      4. In the absence of a free aldehyde, the molecule must be cyclic. However, if there is an ester 
         (lactone) fingerprint—indicated by a pattern "[CX3](=O)[O]"—the structure is rejected.
         (Cyclic aldoses should display a hemiacetal functionality rather than a lactone.)
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule is an aldopentose; False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # 1. Verify that the longest carbon chain among carbon atoms is 5.
    chain_length = longest_carbon_chain(mol)
    if chain_length != 5:
        return False, f"Longest contiguous carbon chain is {chain_length}; expected 5 for an aldopentose"
    
    # 2. Count oxygen atoms.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 5:
        return False, f"Number of oxygen atoms is {oxygen_count}; expected 5 for C5H10O5"
    
    # 3. Define SMARTS patterns.
    # Aldehyde pattern: carbonyl carbon with one hydrogen.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    # Lactone (ester) pattern: carbonyl carbon connected to an oxygen (without an H).
    lactone_pattern = Chem.MolFromSmarts("[CX3](=O)[O]")
    
    # 4. Classification based on presence/absence of aldehyde.
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Open-chain aldopentose: aldehyde group detected."
    else:
        # Molecule is cyclic. Many aldopentoses exist in a hemiacetal (cyclic) form.
        if mol.GetRingInfo().NumRings() > 0:
            # If a lactone functionality is detected then likely this is not a true sugar.
            if mol.HasSubstructMatch(lactone_pattern):
                return False, "Cyclic molecule contains lactone functionality; not an aldopentose."
            else:
                return True, "Cyclized aldopentose: potential open-chain aldehyde form upon ring opening."
        else:
            return False, "No aldehyde group detected and molecule is acyclic; does not meet aldopentose criteria."

# (Optional debug/testing calls could go here, but per instructions we provide a single code block.)