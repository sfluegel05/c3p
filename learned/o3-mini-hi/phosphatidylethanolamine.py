"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: Phosphatidylethanolamine (PE)
Definition: A class of glycerophospholipids in which a phosphatidyl group
            is esterified to the hydroxy group of ethanolamine.
Heuristics used (improved):
  1. Identify a PE head group by traversing from a phosphorus atom: look for a path of atoms
     P – O – C – C – N. The nitrogen must either be neutral or carry a +1 formal charge,
     provided it has at least one attached hydrogen (to rule out quaternary ammonium species).
  2. Require the presence of at least one ester bond, defined as the fragment [CX3](=O)[OX2].
     If two or more ester bonds are present then a higher molecular weight threshold (≥400 Da)
     is used to help capture diacyl-PE versus lyso-PE (which may be lower, ≥300 Da).
  3. Verify that at least one phosphorus atom is present.
Note: This heuristic approach is not perfect – molecules can be mis‐classified. The improvements
here allow protonated (but not quaternary) nitrogen in the head group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine (PE) based on its SMILES string.
    Heuristics:
      - The molecule must contain at least one phosphorus atom.
      - There must be a head group fragment following the pattern:
            P – O – C – C – N
        where the nitrogen must either be neutral or if positively charged must have at least one
        attached hydrogen (to avoid quaternary ammonium groups as in phosphocholine).
      - At least one ester bond (SMARTS "[CX3](=O)[OX2]") must be found.
        If two or more ester bonds are present, a higher molecular weight (≥400 Da) is required,
        otherwise a minimal weight of 300 Da is allowed (to capture lyso-PE).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a PE, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Look for the phosphoethanolamine head group.
    # We search for a phosphorus atom (P) connected to an oxygen (O),
    # which is connected to a carbon (c1), then to another carbon (c2), then to a nitrogen (N).
    headgroup_found = False
    for p in mol.GetAtoms():
        if p.GetAtomicNum() != 15:  # Must be phosphorus
            continue
        # Iterate through neighbors: need an oxygen attached to the phosphorus.
        for o in p.GetNeighbors():
            if o.GetAtomicNum() != 8:
                continue
            # Look for a carbon attached to this oxygen (c1)
            for c1 in o.GetNeighbors():
                if c1.GetIdx() == p.GetIdx() or c1.GetAtomicNum() != 6:
                    continue
                # Instead of imposing too strict a neighbor-count check, simply continue the chain.
                for c2 in c1.GetNeighbors():
                    if c2.GetIdx() == o.GetIdx() or c2.GetAtomicNum() != 6:
                        continue
                    # Now from c2, look for a nitrogen atom (N) that is not the same as c1.
                    for n in c2.GetNeighbors():
                        if n.GetIdx() == c2.GetIdx() or n.GetAtomicNum() != 7:
                            continue
                        # Accept the nitrogen if it is either neutral or if it is positively charged
                        # (formal charge +1) but has at least one hydrogen (to avoid quaternary ammonium).
                        n_charge = n.GetFormalCharge()
                        n_hcount = n.GetTotalNumHs()
                        if n_charge == 0 or (n_charge == 1 and n_hcount > 0):
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
        return False, "Missing or incorrect phosphoethanolamine head group pattern (expected P-O-CH2-CH2-N with neutral or protonated N bearing Hs)"
    
    # 2. Check for ester groups.
    # An ester bond is defined as [CX3](=O)[OX2]. We use a SMARTS pattern search.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if ester_pattern is None:
        return False, "Error parsing ester SMARTS pattern"
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    n_esters = len(ester_matches)
    if n_esters < 1:
        return False, f"Insufficient ester groups detected (found {n_esters}; require at least 1)"

    # 3. Verify that at least one phosphorus atom is present.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphorus atom found in structure"

    # 4. Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # If there are two or more ester bonds, then require a higher molecular weight threshold.
    if n_esters >= 2 and mol_wt < 400:
        return False, f"Molecular weight too low for a diacyl phosphatidylethanolamine ({mol_wt:.1f} Da; expected >= 400 Da)"
    elif n_esters == 1 and mol_wt < 300:
        return False, f"Molecular weight too low for a lyso-phosphatidylethanolamine ({mol_wt:.1f} Da; expected >= 300 Da)"
    
    return True, f"Structure contains a phosphoethanolamine head group (P-O-CH2-CH2-N), {n_esters} ester bond(s), phosphorus, and acceptable molecular weight ({mol_wt:.1f} Da)"

# Simple test routine if run as a script (can be removed or commented when used as a module)
if __name__ == "__main__":
    # Example SMILES from one of the given examples:
    test_smiles = "P(OCC(OC(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O"
    result, reason = is_phosphatidylethanolamine(test_smiles)
    print(result, reason)