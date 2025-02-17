"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: Sulfolipid
A compound containing a sulfonic acid residue joined (directly or via a single oxygen) to a long lipid chain.
The algorithm checks:
  - The molecule must be above 300 Da.
  - It must exhibit a sulfonic acid group (S with two double-bonded oxygens and one single-bonded O).
  - One of the non-oxygen substituents (or an O-bridged neighbor) of that S must be a sp3 carbon
    (not in a ring) that supports a long, uninterrupted aliphatic chain (≥12 carbons) when traversed by DFS.
If these conditions are met, the molecule is classified as a sulfolipid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is defined as a compound containing a sulfonic acid residue attached (directly or via a single oxygen)
    to a long aliphatic lipid chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a sulfolipid, False otherwise.
        str: Reason for the classification result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low for a sulfolipid (mol wt = {mol_wt:.1f})"
    
    # A helper routine to compute the length of a continuous non-cyclic, aliphatic chain.
    def dfs_chain(atom, visited):
        count = 1
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in visited:
                continue
            # Only traverse sp3 carbons that are non-aromatic and not in any ring.
            if (nbr.GetAtomicNum() == 6 and 
                nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and
                not nbr.GetIsAromatic() and
                not nbr.IsInRing()):
                new_visited = visited.copy()
                new_visited.add(nbr.GetIdx())
                candidate = 1 + dfs_chain(nbr, new_visited)
                if candidate > count:
                    count = candidate
        return count

    # Define SMARTS patterns for sulfonic acid – it may be depicted as S(=O)(=O)[O-] or S(=O)(=O)O.
    sulfonate_smarts_1 = Chem.MolFromSmarts("S(=O)(=O)[O-]")
    sulfonate_smarts_2 = Chem.MolFromSmarts("S(=O)(=O)O")
    
    sulfonate_atoms = []  # will hold sulfur atoms that match sulfonate patterns.
    for smarts in [sulfonate_smarts_1, sulfonate_smarts_2]:
        matches = mol.GetSubstructMatches(smarts, useChirality=True)
        for match in matches:
            # match is a tuple of atom indices; for our patterns the first atom is the sulfur.
            sulfonate_atoms.append(match[0])
    # Remove duplicates
    sulfonate_atoms = set(sulfonate_atoms)
    
    if not sulfonate_atoms:
        return False, "No sulfonic acid group detected"
    
    # For each sulfonate sulfur atom, check its connectivity for a lipid chain.
    for s_idx in sulfonate_atoms:
        atom_S = mol.GetAtomWithIdx(s_idx)
        # We search among the neighbors for a candidate carbon.
        candidate_c_atoms = []
        for nbr in atom_S.GetNeighbors():
            # if directly attached carbon and not oxygen
            if nbr.GetAtomicNum() == 6:
                candidate_c_atoms.append(nbr)
            # if neighbor is oxygen (with a single bond), follow that oxygen to see if it connects to a carbon.
            elif nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom_S.GetIdx(), nbr.GetIdx())
                if bond and abs(bond.GetBondTypeAsDouble() - 1.0) < 0.01:
                    for nn in nbr.GetNeighbors():
                        if nn.GetAtomicNum() == 6 and nn.GetIdx() != nbr.GetIdx():
                            candidate_c_atoms.append(nn)
        # Evaluate each candidate carbon:
        for c_atom in candidate_c_atoms:
            # Exclude candidate carbons that are in rings or not sp3.
            if c_atom.IsInRing() or c_atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            chain_length = dfs_chain(c_atom, {c_atom.GetIdx()})
            if chain_length >= 12:
                reason = (f"Found sulfonic acid group attached to an aliphatic chain (chain length = {chain_length}) "
                          f"via carbon atom idx {c_atom.GetIdx()}")
                return True, reason
            # Else record if chain is found but too short
            short_reason = f"Candidate chain attached has only {chain_length} carbons"
        # If we get here for this sulfur, none of its candidate carbons yielded a long chain.
    return False, "No sulfonic acid residue found attached to a sufficiently long aliphatic chain"

# Example usage (for testing):
if __name__ == '__main__':
    test_smiles = [
        "CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H](N)CO[C@@H]1O[C@H](COS(O)(=O)=O)[C@H](O)[C@H](O)[C@H]1O",  # psychosine sulfate
        "CC1=CC=C(C=C1)S(=O)(=O)O",  # a generic aromatic sulfonate (should not be sulfolipid)
    ]
    for smi in test_smiles:
        is_sulfo, reason = is_sulfolipid(smi)
        print(f"SMILES: {smi}\nClassified as sulfolipid? {is_sulfo}\nReason: {reason}\n")