"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: Sulfolipid
A compound containing a sulfonic acid residue (S(=O)(=O)[O-] or S(=O)(=O)O) attached (directly or via a single oxygen)
to a long, uninterrupted aliphatic chain. The algorithm checks:
  - Molecular weight must be above 300 Da.
  - The molecule must contain a sulfonic acid group.
  - One of the substituents (or one linked via one oxygen) of that sulfur must be a sp³ carbon (non‐aromatic, not in a ring)
    that supports a continuous chain of at least 12 sp³ carbons.
If these conditions are met, the molecule is classified as a sulfolipid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is defined as a compound containing a sulfonic acid residue attached (directly or via a single oxygen)
    to a long, uninterrupted aliphatic lipid chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a sulfolipid, False otherwise.
        str: Reason for the classification result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that the molecular weight is high enough for a lipid
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low for a sulfolipid (mol wt = {mol_wt:.1f})"

    # Helper: a DFS to find the maximum length of continuous aliphatic (sp3, non-aromatic, acyclic) carbon chain.
    def dfs_chain_length(atom, visited):
        max_length = 1  # count the current carbon
        for nbr in atom.GetNeighbors():
            # Only traverse through carbon atoms meeting our aliphatic criteria.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                if (nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and 
                    not nbr.GetIsAromatic() and 
                    not nbr.IsInRing()):
                    new_visited = visited.copy()
                    new_visited.add(nbr.GetIdx())
                    length = 1 + dfs_chain_length(nbr, new_visited)
                    if length > max_length:
                        max_length = length
        return max_length

    # Define SMARTS for the sulfonic acid group.
    sulfonate_smarts_1 = Chem.MolFromSmarts("S(=O)(=O)[O-]")  # negatively charged sulfonate
    sulfonate_smarts_2 = Chem.MolFromSmarts("S(=O)(=O)O")       # neutral sulfonic acid

    sulfonate_atoms = set()
    for smarts in [sulfonate_smarts_1, sulfonate_smarts_2]:
        matches = mol.GetSubstructMatches(smarts, useChirality=True)
        for match in matches:
            sulfonate_atoms.add(match[0])  # The sulfur is the first atom in the pattern

    if not sulfonate_atoms:
        return False, "No sulfonic acid group detected"

    # Now for each sulfonate sulfur, check each neighbor for a candidate aliphatic chain.
    for s_idx in sulfonate_atoms:
        sulfur_atom = mol.GetAtomWithIdx(s_idx)
        candidate_c_atoms = []
        # Look at atoms directly bonded to sulfur
        for nbr in sulfur_atom.GetNeighbors():
            # Case 1: directly attached carbon
            if nbr.GetAtomicNum() == 6:
                candidate_c_atoms.append(nbr)
            # Case 2: when an oxygen bridges between sulfur and a carbon
            elif nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(sulfur_atom.GetIdx(), nbr.GetIdx())
                if bond and abs(bond.GetBondTypeAsDouble() - 1.0) < 0.01:
                    # For each neighbor of this oxygen (except the sulfur), check if it's a carbon.
                    for nn in nbr.GetNeighbors():
                        if nn.GetAtomicNum() == 6 and nn.GetIdx() != sulfur_atom.GetIdx():
                            candidate_c_atoms.append(nn)

        # Evaluate candidates: they must be sp3, non-aromatic, and not in a ring.
        for cand in candidate_c_atoms:
            if cand.IsInRing() or cand.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue  # skip unsuitable candidates
            # Compute the length of the continuous carbon chain from this candidate.
            chain_length = dfs_chain_length(cand, {cand.GetIdx()})
            if chain_length >= 12:
                reason = (f"Found sulfonic acid group attached to an aliphatic chain (chain length = {chain_length}) "
                          f"via carbon atom idx {cand.GetIdx()}")
                return True, reason
    # If no candidate carbon attached (or via a single oxygen) led to a long chain, return False.
    return False, "No sulfonic acid residue found attached to a sufficiently long aliphatic chain"

# Example usage if running this script standalone (for testing):
if __name__ == '__main__':
    test_smiles = [
        # psychosine sulfate
        "CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H](N)CO[C@@H]1O[C@H](COS(O)(=O)=O)[C@H](O)[C@H](O)[C@H]1O",
        # A generic aromatic sulfonate (not a sulfolipid)
        "CC1=CC=C(C=C1)S(=O)(=O)O",
    ]
    for smi in test_smiles:
        is_sulfo, reason = is_sulfolipid(smi)
        print(f"SMILES: {smi}\nClassified as sulfolipid? {is_sulfo}\nReason: {reason}\n")