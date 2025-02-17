"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
Definition: An alkylamine in which one or more alkyl substituents bears an aromatic group.
Operationally, we require that at least one non‐aromatic amine nitrogen (i.e. not in any aromatic ring)
has at least one carbon-only substituent that (via a chain of single C–C bonds, all using sp3 carbons)
leads to an aromatic carbon. To avoid spurious results we ignore any direct (distance 1) attachment
of an aromatic carbon (which would be an arylamine) and only accept a branching if an aromatic carbon
appears at a bond distance of at least 2 (and at most 3 bonds away from the nitrogen).
"""

from rdkit import Chem
from collections import deque

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    
    Approach:
      1. Parse the SMILES string.
      2. For each nitrogen (atomic number 7) that is not aromatic, examine each substituent.
      3. Only consider substituents that are attached by a single bond and begin with a carbon.
         If the immediate carbon is aromatic, skip it (direct attachment would be an arylamine).
         Also require that that starting carbon is sp3 (as required for an alkyl chain).
      4. Use breadth‐first search (BFS) along the branch – but only traverse through carbon atoms that are:
            (a) non‐aromatic, and
            (b) sp3 hybridized.
         (Aromatic carbons are only “accepted” if they are encountered as a terminal hit and are at least 2 bonds 
         away from the nitrogen.)
      5. The search is capped at a maximum depth (here 3 bonds including the first C–N bond).
      6. If any branch meets the criteria we return True along with a reason; otherwise we return False.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): Tuple with boolean classification and a message explaining the decision.
      If the SMILES cannot be parsed, returns (False, "Invalid SMILES string").
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    max_depth = 3  # maximum total bond distance allowed from the N atom
    
    # Iterate over all atoms looking for non‐aromatic amine nitrogen atoms.
    for n_atom in mol.GetAtoms():
        if n_atom.GetAtomicNum() != 7:
            continue
        if n_atom.GetIsAromatic():
            continue  # skip if the N is aromatic
        n_idx = n_atom.GetIdx()
        
        # Examine each bond of the nitrogen:
        for bond in n_atom.GetBonds():
            # Only consider single bonds (alkyl bonds)
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            neighbor = bond.GetOtherAtom(n_atom)
            # We only consider branches that “start” with a carbon atom.
            if neighbor.GetAtomicNum() != 6:
                continue
            # If the atom immediately attached is aromatic, that is a direct N–aryl bond;
            # skip this branch because we require the aromatic substituent to be on an alkyl branch.
            if neighbor.GetIsAromatic():
                continue
            # Also require that the branch carbon is sp3 (to ensure an alkyl chain)
            if neighbor.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            
            # Start a BFS from this neighbor along the alkyl substituent.
            # The distance here counts bonds from the nitrogen.
            start_idx = neighbor.GetIdx()
            visited = {n_idx, start_idx}
            # Start at depth 1. (Benzylamine, N-CH2-Ph, would show the aromatic ring at depth 2.)
            queue = deque([(start_idx, 1)])
            
            while queue:
                current_idx, depth = queue.popleft()
                current_atom = mol.GetAtomWithIdx(current_idx)
                # If we are at a distance of at least 2 and current atom is aromatic, then we have found
                # an aralkyl substituent.
                if depth >= 2 and current_atom.GetIsAromatic():
                    reason = (f"Found a non‐aromatic amine nitrogen (atom {n_idx}) with an aromatic substituent "
                              f"(atom {current_idx}) at a bond distance of {depth}. Molecule classified as aralkylamine.")
                    return True, reason
                # Do not extend beyond max_depth
                if depth >= max_depth:
                    continue
                # Expand allowed branches: only traverse single bonds to carbon atoms that are non‐aromatic and sp3.
                for nbr in current_atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in visited:
                        continue
                    bond2 = mol.GetBondBetweenAtoms(current_idx, nbr_idx)
                    if bond2.GetBondType() != Chem.rdchem.BondType.SINGLE:
                        continue
                    # Only allow the branch to continue if the neighbor is carbon.
                    if nbr.GetAtomicNum() != 6:
                        continue
                    new_depth = depth + 1
                    # If neighbor is aromatic and is reached at a valid distance, we accept it.
                    if nbr.GetIsAromatic():
                        if new_depth >= 2 and new_depth <= max_depth:
                            reason = (f"Found a non‐aromatic amine nitrogen (atom {n_idx}) with an aromatic substituent "
                                      f"(atom {nbr_idx}) at a bond distance of {new_depth}. Molecule classified as aralkylamine.")
                            return True, reason
                        continue
                    # For non‐aromatic carbons, we require the carbon to be sp3 to consider it part of an alkyl chain.
                    if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                        continue
                    # Only traverse further if we are within the allowed depth.
                    if new_depth < max_depth:
                        visited.add(nbr_idx)
                        queue.append((nbr_idx, new_depth))
                    # If new_depth equals max_depth we do not extend further (but we already checked for aromatic).
            # End BFS for this branch.
    # If no branch meets the criteria, classify as not an aralkylamine.
    return (False, 
            "No aralkylamine substructure found: "
            "no non‐aromatic amine nitrogen is connected via a suitable alkyl chain "
            "to an aromatic carbon within allowed bond distance.")

# Optional test block; remove or comment this block out when deploying as a module.
if __name__ == "__main__":
    # Here are some example test cases, including several true and false positives/negatives.
    test_cases = [
        ("NCc1ccccc1", "benzylamine"),  # expected True: N-CH2 (alkyl) then phenyl at distance 2.
        ("NCCc1ccccc1", "phenethylamine"),  # expected True: chain length 2 reaching an aromatic ring.
        ("c1ccc(N)cc1", "aniline (should not classify, N aromatic)"),  # expected False.
        ("OCCNC1=CC=CC=C1", "2-Anilinoethanol"),  # expected True: branch from N to CH2 then aromatic.
        ("Cl.C1=CC=CC(=C1)C(C2CCCCC2)(CCN3CCCCC3)O", "Trihexyphenidyl hydrochloride"),  # expected True.
    ]
    
    for smi, name in test_cases:
        res, explanation = is_aralkylamine(smi)
        print(f"SMILES: {smi}\nMolecule: {name}\nClassification: {res}\nReason: {explanation}\n")