"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: Phosphatidylethanolamine
Definition: A class of glycerophospholipids in which a phosphatidyl group is esterified 
to the hydroxy group of ethanolamine.
Heuristics used (improved):
  - Instead of relying solely on a fixed SMARTS, in the head‐group search we traverse from a
    phosphorus atom. For each phosphorus we examine its oxygen neighbors and then search for a two‐carbon
    chain ending in a nitrogen. That nitrogen must be neutral (formal charge 0) so that quaternary ammoniums
    (e.g. phosphocholine) are excluded.
  - We then require at least one ester bond – defined as a “[CX3](=O)[OX2]” fragment –
    though if two are present we require a higher molecular weight (≥ 400 Da) than if only one is found (≥ 300 Da).
  - In addition, we verify that at least one phosphorus atom is present.
Note: This heuristic approach is not perfect, but we hope that by relaxing the rules for the acyl chains
(a bit to capture lyso-PE) and checking the formal charge of the amine the overall classification F1 will improve.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine (PE) based on its SMILES string.
    A PE is characterized by:
      1. A phospho head group esterified to an ethanolamine-like fragment (specifically,
         a P–O–CH2–CH2–N fragment; if the nitrogen is substituted it must be neutral, to exclude PC).
      2. At least one ester bond ([CX3](=O)[OX2]) is present.
         If two ester bonds are present, a higher molecular weight (≥ 400 Da) is required.
         (This is to try to capture lyso-PE cases with only one fatty acyl chain.)
      3. The molecule must contain at least one phosphorus atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a phosphatidylethanolamine, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Look for a phosphoethanolamine head group. We iterate over phosphorus atoms.
    headgroup_found = False
    for p_atom in mol.GetAtoms():
        if p_atom.GetAtomicNum() != 15:  # phosphorus
            continue
        # For each oxygen neighbor attached to phosphorus
        for o_atom in p_atom.GetNeighbors():
            if o_atom.GetAtomicNum() != 8:
                continue
            # Look for a carbon (first carbon in the chain) attached to the oxygen (but not back to phosphorus)
            for c1 in o_atom.GetNeighbors():
                if c1.GetIdx() == p_atom.GetIdx() or c1.GetAtomicNum() != 6:
                    continue
                # We expect c1 (should be like CH2) to have 2 heavy-atom neighbors in an ideal ETHANOLAMINE fragment.
                heavy_neighbors_c1 = [a for a in c1.GetNeighbors() if a.GetAtomicNum() > 1]
                if len(heavy_neighbors_c1) != 2:
                    continue
                # Now from c1, search for a second carbon (c2) that is not the oxygen.
                for c2 in c1.GetNeighbors():
                    if c2.GetIdx() == o_atom.GetIdx() or c2.GetAtomicNum() != 6:
                        continue
                    heavy_neighbors_c2 = [a for a in c2.GetNeighbors() if a.GetAtomicNum() > 1]
                    if len(heavy_neighbors_c2) != 2:
                        continue
                    # Look for a nitrogen neighbor on c2 (other than c1)
                    for n_atom in c2.GetNeighbors():
                        if n_atom.GetIdx() == c1.GetIdx():
                            continue
                        if n_atom.GetAtomicNum() != 7:
                            continue
                        # Check that the nitrogen is neutral (formal charge 0)
                        if n_atom.GetFormalCharge() != 0:
                            continue
                        headgroup_found = True
                        break
                    if headgroup_found:
                        break
                if headgroup_found:
                    break
            if headgroup_found:
                break
        if headgroup_found:
            break
                    
    if not headgroup_found:
        return False, "Missing or incorrect phosphoethanolamine head group pattern (expected P–O–CH2–CH2–N with neutral N)"
    
    # 2. Check for ester groups.
    # Define an ester bond as [CX3](=O)[OX2] (i.e. a carbonyl carbon attached to an oxygen).
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if ester_pattern is None:
        return False, "Error parsing the ester SMARTS pattern"
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    n_esters = len(ester_matches)
    
    if n_esters < 1:
        return False, f"Insufficient ester groups detected (found {n_esters}; require at least 1)"
    
    # 3. Verify the presence of at least one phosphorus atom.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphorus atom found in the structure"
    
    # 4. Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Use different thresholds: if two or more esters, require >= 400 Da; if one ester (lyso-PE) allow lower threshold (>=300 Da)
    if n_esters >= 2 and mol_wt < 400:
        return False, f"Molecular weight too low for a diacyl phosphatidylethanolamine ({mol_wt:.1f} Da; expected >= 400 Da)"
    elif n_esters == 1 and mol_wt < 300:
        return False, f"Molecular weight too low for a lyso-phosphatidylethanolamine ({mol_wt:.1f} Da; expected >= 300 Da)"
    
    return True, f"Structure contains a phosphoethanolamine head group (P-O-CH2-CH2-N with neutral N), {n_esters} ester bond(s), phosphorus, and acceptable molecular weight ({mol_wt:.1f} Da)"

# Simple test routine if run as a script (this section can be removed or commented out when used as a module)
if __name__ == "__main__":
    # Example SMILES taken from one of the provided examples.
    test_smiles = "P(OCC(OC(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O"
    result, reason = is_phosphatidylethanolamine(test_smiles)
    print(result, reason)