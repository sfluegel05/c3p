"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: Lipid Hydroperoxide – any lipid carrying one or more hydroperoxy (-OOH) substituents.
Improved heuristics:
  1. Look for a neutral hydroperoxy group using a SMARTS pattern that matches any non‐aromatic carbon attached to an O–O–H group.
  2. Check that the atoms in the found group and their immediate neighbors are all formally neutral.
  3. Evaluate that the hydroperoxy substituent is attached to a long contiguous aliphatic chain (>=6 carbons) 
     by traversing non‐aromatic, non‐ring carbon atoms.
  4. Ensure the overall molecule is uncharged (formal charge 0) and has a reasonable molecular weight (>250 Da).
If all conditions are met, the molecule is classified as a lipid hydroperoxide.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule represents a lipid hydroperoxide.
    Here a lipid hydroperoxide is defined as a molecule that carries at least one neutral hydroperoxy (-OOH)
    substituent attached to a contiguous (>= 6 carbons) aliphatic (non-aromatic, non-ring) segment,
    with the molecule overall uncharged and relatively large (MW > 250 Da).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a lipid hydroperoxide; False otherwise.
        str: Explanation for the classification decision.
    """
    # 1. Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1a. Reject if the overall molecule has a non-zero net formal charge.
    if Chem.GetFormalCharge(mol) != 0:
        return False, f"Molecule has a net formal charge of {Chem.GetFormalCharge(mol)}, expected 0"
    
    # 2. Look for a hydroperoxy group.
    # We use a looser SMARTS: any non-aromatic carbon ([C;!a]) attached to an O–O–H group.
    hydroperoxy_pattern = Chem.MolFromSmarts("[C;!a]-[OX2]-[OX2H]")
    hydro_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    valid_hydro_match = None
    
    def local_environment_is_neutral(atom_idx):
        """Check that atom and its immediate neighbors are all formally neutral."""
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetFormalCharge() != 0:
            return False
        for nb in atom.GetNeighbors():
            if nb.GetFormalCharge() != 0:
                return False
        return True

    # Look for a matching hydroperoxy group whose atoms and their neighbors are all neutral.
    for match in hydro_matches:
        # match: (hydrocarbon, middle oxygen, terminal oxygen with H)
        if (local_environment_is_neutral(match[0]) and 
            local_environment_is_neutral(match[1]) and 
            local_environment_is_neutral(match[2])):
            valid_hydro_match = match
            break
    if valid_hydro_match is None:
        return False, "No valid neutral hydroperoxy (-OOH) group found"
    
    # 3. Check for a contiguous aliphatic chain attached to the hydroperoxy carbon.
    # Starting from the hydroperoxy carbon (the first atom in the match), perform a DFS.
    # We only traverse carbon atoms that are non-aromatic and not part of any ring.
    def longest_chain_depth(atom_idx, visited):
        """Recursively compute the maximum contiguous chain length from a given carbon atom.
           Only counts neighbors that are carbons, non-aromatic, and not in a ring.
           Returns the count including the starting atom.
        """
        visited.add(atom_idx)
        current_atom = mol.GetAtomWithIdx(atom_idx)
        max_depth = 1  # count the current atom
        for neighbor in current_atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if (neighbor.GetAtomicNum() == 6 and
                not neighbor.GetIsAromatic() and
                not neighbor.IsInRing() and
                n_idx not in visited):
                depth = 1 + longest_chain_depth(n_idx, visited.copy())
                if depth > max_depth:
                    max_depth = depth
        return max_depth

    hydro_carbon_idx = valid_hydro_match[0]
    chain_length = longest_chain_depth(hydro_carbon_idx, set())
    if chain_length < 6:
        return False, f"Longest contiguous aliphatic chain (length {chain_length}) attached to the hydroperoxy group is too short"

    # 4. Check overall molecular weight; many lipids have MW > 250 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical lipid"
    
    return True, "Molecule is classified as a lipid hydroperoxide: contains a neutral hydroperoxy (-OOH) group attached to a long aliphatic chain"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples (5S-HpEPE).
    test_smiles = "O(O)[C@@H](CCCC(O)=O)/C=C/C=C\\C/C=C\\C/C=C\\CC"
    result, reason = is_lipid_hydroperoxide(test_smiles)
    print(result, reason)