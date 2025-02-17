"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: Lipid Hydroperoxide – any lipid carrying one or more hydroperoxy (-OOH) substituents.
Heuristics:
  1. Look for at least one hydroperoxy group using a SMARTS pattern that requires the pattern [#6]-[OX2]-[OX2H]
     (i.e. a carbon directly bonded to an oxygen, which is then bonded to an oxygen bearing a hydrogen).
     Also verify that these atoms are formally neutral to avoid matching ionized groups.
  2. From the carbon that bears the hydroperoxy group, check if at least one contiguous chain of non‐aromatic,
     non‐ring carbon atoms (using a DFS) has a length of at least 8.
  3. Check that the overall molecular weight is above a modest cutoff.
If all conditions are met, the molecule is classified as a lipid hydroperoxide.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule represents a lipid hydroperoxide.
    A lipid hydroperoxide is defined here as a molecule that contains a hydroperoxy (-OOH)
    substituent attached to a long, contiguous unbranched aliphatic chain (>=8 carbons),
    and which has a molecular weight above a minimal threshold.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a lipid hydroperoxide; False otherwise.
        str: Explanation of the classification decision.
    """
    # 1. Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Look for a hydroperoxy group.
    # Here we require a pattern: carbon attached to an oxygen, connected to an oxygen with hydrogen.
    # SMARTS: [#6]-[OX2]-[OX2H]
    hydroperoxy_pattern = Chem.MolFromSmarts("[#6]-[OX2]-[OX2H]")
    hydro_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    valid_hydro_match = None
    for match in hydro_matches:
        # match is a tuple of atom indices: (carbon, oxygen, terminal oxygen (with H))
        # Ensure that each of these atoms has formal charge 0, so that deprotonated groups do not qualify.
        if all(mol.GetAtomWithIdx(idx).GetFormalCharge() == 0 for idx in match):
            valid_hydro_match = match
            break
    if valid_hydro_match is None:
        return False, "No valid (neutral) hydroperoxy (-OOH) group found"
    
    # 3. Check that there is a long aliphatic chain.
    # We want to ensure that from the carbon (first atom in the hydroperoxy match)
    # one can traverse a chain of contiguous carbon atoms (non-aromatic, non-ring)
    # of length at least 8.
    def longest_chain_depth(atom_idx, visited):
        """Recursively determine the maximum chain length from a starting carbon atom.
           Only non-aromatic carbons that are not in a ring are traversed.
           Returns the chain length including the starting atom.
        """
        visited.add(atom_idx)
        current_atom = mol.GetAtomWithIdx(atom_idx)
        max_depth = 1  # count current atom
        # Iterate over neighbors that are carbon atoms.
        for neighbor in current_atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            # conditions: must be carbon, not visited already, non-aromatic and not in a ring.
            if neighbor.GetAtomicNum() == 6 and (not neighbor.GetIsAromatic()) and (not neighbor.IsInRing()) and (n_idx not in visited):
                depth = 1 + longest_chain_depth(n_idx, visited.copy())
                if depth > max_depth:
                    max_depth = depth
        return max_depth

    # Use the carbon atom from the valid hydroperoxy match as a starting point.
    hydro_carbon_idx = valid_hydro_match[0]
    chain_length = longest_chain_depth(hydro_carbon_idx, set())
    if chain_length < 8:
        return False, f"Longest contiguous aliphatic chain (length {chain_length}) attached to hydroperoxy group is too short"
    
    # 4. Check overall molecular weight; many lipids have MW > 250 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for typical lipid"
    
    return True, "Molecule is classified as a lipid hydroperoxide: contains a neutral hydroperoxy (-OOH) group attached to a long aliphatic chain"

# Example usage:
if __name__ == "__main__":
    # Testing one of the known true positives, for example 5S-HpEPE:
    test_smiles = "O(O)[C@@H](CCCC(O)=O)/C=C/C=C\\C/C=C\\C/C=C\\CC"
    result, reason = is_lipid_hydroperoxide(test_smiles)
    print(result, reason)