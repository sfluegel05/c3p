"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: Phosphatidylethanolamine
Definition: A class of glycerophospholipids in which a phosphatidyl group is esterified 
to the hydroxy group of ethanolamine.
Heuristics used (improved):
  - Instead of a fixed SMARTS, we manually walk from a phosphorus atom. For each phosphorus,
    we check if one of its oxygen neighbors is bound to a two‐carbon chain that terminates
    in a nitrogen. In this chain, the first carbon (directly connected to the oxygen) must have
    exactly two heavy‐atom neighbors and the second carbon (connected to the first and the nitrogen)
    must also have exactly two heavy‐atom neighbors.  This tends to capture the ethanolamine 
    head group – CH2–CH2–NH2 (or N-alkylated variants) – while rejecting similar motifs (for example,
    those in phosphatidylserine which have an extra substituent).
  - Require at least two ester bonds (each detected as [CX3](=O)[OX2]).
  - Verify that a phosphorus atom is present.
  - Enforce a molecular weight threshold (≥ 400 Da).
Note: This heuristic approach is not perfect, but should improve on our previous attempt.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine (PE) based on its SMILES string.
    A PE is characterized by:
      1. A phospho head group esterified to an ethanolamine-like fragment (specifically,
         a P–O–CH2–CH2–N fragment where the CH2 groups are unadorned, ensuring it is not serine).
      2. At least two ester bonds (indicating two fatty acyl chains).
      3. The presence of at least one phosphorus atom.
      4. A molecular weight typically above 400 Da.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a phosphatidylethanolamine, False otherwise.
        str: Explanation for the classification result.
    """
    
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the phosphoethanolamine head group.
    # Instead of relying solely on SMARTS, we search for a phosphorus atom (P)
    # that is linked to an oxygen, which in turn bonds to a two‐carbon chain ending in a nitrogen.
    headgroup_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # Phosphorus
            # For each neighbor O attached to the phosphorus
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() != 8:
                    continue
                # Now from this oxygen, get its other neighbor(s) (excluding our P)
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() == atom.GetIdx():
                        continue
                    if nbr2.GetAtomicNum() != 6:  # Must be carbon (first carbon in chain)
                        continue
                    # Check that this first carbon (c1) appears to be a CH2-like fragment:
                    # we require that it has exactly two heavy-atom neighbors (typically the oxygen and the next carbon)
                    if len([a for a in nbr2.GetNeighbors() if a.GetAtomicNum() > 1]) != 2:
                        continue
                    # Now look for the second carbon attached to c1 (other than the oxygen)
                    for nbr3 in nbr2.GetNeighbors():
                        if nbr3.GetIdx() == nbr.GetIdx():
                            continue
                        if nbr3.GetAtomicNum() != 6:
                            continue
                        # Check that the second carbon (c2) is only connected to c1 and one more heavy atom (the nitrogen)
                        if len([a for a in nbr3.GetNeighbors() if a.GetAtomicNum() > 1]) != 2:
                            continue
                        # Finally, check for an adjacent nitrogen (this can be ammonium or substituted amine)
                        for nbr4 in nbr3.GetNeighbors():
                            if nbr4.GetIdx() == nbr2.GetIdx():
                                continue
                            if nbr4.GetAtomicNum() == 7:  # nitrogen found
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
        return False, "Missing or incorrect phosphoethanolamine head group pattern"
    
    # 2. Check for at least two ester bonds.
    # An ester bond is defined as [CX3](=O)[OX2].
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if ester_pattern is None:
        return False, "Error parsing ester SMARTS pattern"
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester groups for fatty acyl chains (found {len(ester_matches)}, require at least 2)"
    
    # 3. Verify the presence of at least one phosphorus atom.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphorus atom found in the structure"
    
    # 4. Check molecular weight – PE molecules are typically heavy.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low for a phosphatidylethanolamine ({mol_wt:.1f} Da; expected >= 400 Da)"
    
    return True, "Structure contains a phosphoethanolamine head group (P-O-CH2-CH2-N), at least 2 ester bonds, phosphorus, and acceptable molecular weight"

# For testing purposes (this section can be removed or commented out when used as a module)
if __name__ == "__main__":
    # Example test SMILES (taken from one of the provided examples)
    test_smiles = "P(OCC(OC(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O"
    result, reason = is_phosphatidylethanolamine(test_smiles)
    print(result, reason)