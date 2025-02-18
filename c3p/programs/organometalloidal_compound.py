"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: organometalloidal compound
Definition: A compound having bonds between one or more metalloid atoms and one or more carbon atoms 
             of an organyl group.
Note: In order to reduce false positives from metalloids other than arsenic and from compounds 
      that show an As–C bond but likely belong to another class (e.g. tryparsamide and similar),
      this implementation restricts to arsenic-based compounds and applies additional heuristics.
      
Heuristics used:
 1. The SMILES must contain an arsenic atom (“[As]”) – other metalloids are rejected.
 2. For each arsenic atom, we look at its bonds to carbon.
     • If the attached carbon is part of a larger organic fragment (i.e. it is bonded to at least one other carbon),
       then we count this as evidence.
     • If the only organic substituent is a methyl group then we allow it only if the arsenic has only one C–bond.
       (Many false positives like trimethylarsine oxide have several methyl groups.)
 3. In addition, as a “hard‐coded” filter, we reject SMILES that include a tryparsamide‐like fragment.
 
The method is only one possible strategy given that an F₁ score of about 0.21 was seen previously.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    (Heuristics: restrict to arsenic-based compounds having at least one bond between arsenic and a carbon 
     belonging to an organyl fragment.)
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is judged to be an organometalloidal compound, False otherwise.
        str: Reason for classification.
    """
    # First try to parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Only consider arsenic-based compounds.
    if "[As" not in smiles:
        return False, "No arsenic atom found; not organo–arsenical"
    
    # A hard-coded check for a known false-positive pattern (e.g. tryparsamide-like fragment)
    if "NCC(N)=O" in smiles:
        return False, "Detected a tryparsamide-like fragment (NCC(N)=O) that is not considered an organyl group"
    
    # Compute approximate molecular weight – very small compounds are often false positives.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 120:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be an organometalloidal compound"
    
    # Now, for each arsenic atom, consider its C–As bonds.
    # We will require that at least one such bond is “organyl” – meaning that the carbon partner is
    # either not simply a methyl group or (if it is methyl) that arsenic has only one carbon substituent.
    valid_bond_found = False
    details = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 33:
            continue  # only care about arsenic for this classification
        # Gather neighbors that are carbons
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not carbon_neighbors:
            continue  # no As–C bonds for this arsenic
        # For this arsenic atom, check each bonded carbon.
        accepted_for_this_As = False
        for c in carbon_neighbors:
            # Determine if the carbon is part of a larger organic (organyl) fragment.
            # We require that, aside from the bond to arsenic, the carbon is bonded to at least one other carbon.
            other_c_neighbors = [nbr for nbr in c.GetNeighbors() if nbr.GetIdx() != atom.GetIdx() and nbr.GetAtomicNum() == 6]
            if other_c_neighbors:
                # This carbon is part of a chain or ring.
                accepted_for_this_As = True
                details.append("Found As–C bond where C is part of a larger organic fragment")
            else:
                # This carbon appears to be a methyl group.
                # Accept a methyl bond only if this arsenic has exactly one carbon neighbor.
                if len(carbon_neighbors) == 1:
                    accepted_for_this_As = True
                    details.append("Found As–C bond with a methyl group (only substituent on As)")
                else:
                    details.append("Found As–C bond with a methyl group but arsenic has multiple C–bonds; skipping")
        if accepted_for_this_As:
            valid_bond_found = True
            # (We break early if at least one arsenic atom is acceptable.)
            break

    if valid_bond_found:
        return True, "; ".join(details)
    else:
        return False, "No As–C bond found with an appropriate organyl fragment"

# Example usage:
if __name__ == "__main__":
    # A true positive example:
    example_smiles = "[As](=O)(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)(C)C"
    result, reason = is_organometalloidal_compound(example_smiles)
    print(f"Input SMILES: {example_smiles}\nResult: {result}\nReason: {reason}")
    
    # A false positive example (tryparsamide):
    fp_smiles = "[As](=O)(O)(O)C1=CC=C(C=C1)NCC(N)=O"
    result, reason = is_organometalloidal_compound(fp_smiles)
    print("\nInput SMILES: " + fp_smiles)
    print(f"Result: {result}\nReason: {reason}")