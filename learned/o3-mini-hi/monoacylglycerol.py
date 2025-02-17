"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: Monoacylglycerol
Definition: A glyceride in which any one of the three hydroxyl groups (on a contiguous, acyclic three–carbon glycerol backbone) is esterified (i.e. bears an –O–C(=O)–R group) while the other two bear free –OH.
Heuristic approach (improved):
  1. Parse the molecule and add explicit hydrogens.
  2. Check that the molecular weight is above a threshold (here 140 Da).
  3. Enumerate all acyclic (non–ring) three–carbon chains (paths of length 2) where each atom is sp³ carbon.
  4. For each such candidate “glycerol” chain, for each carbon, consider the oxygen neighbors not belonging to the chain.
       • Classify an oxygen as “ester” if it is bonded to a carbonyl carbon (i.e. a neighboring carbon has a double bond to oxygen).
       • Classify an oxygen as a free –OH if it carries an explicit hydrogen.
  5. Require that (within one candidate chain) exactly one carbon has an ester oxygen (and no free –OH) and the other two each have at least one free –OH.
  6. If any candidate chain meets these criteria, return True (with explanation); else return False.
  
Note: This is a heuristic method and may not be perfect.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    The method:
      - Adds explicit hydrogens.
      - Checks overall molecular weight (>140 Da).
      - Searches for an acyclic three–carbon chain where one carbon bears an ester (acyl) group and the other two each bear a free hydroxyl.
    
    Args:
        smiles (str): SMILES string.
        
    Returns:
        bool: True if the molecule is classified as a monoacylglycerol.
        str: Reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Work with explicit hydrogens for correct OH detection
    mol = Chem.AddHs(mol)

    # Check molecular weight threshold (e.g., 140 Da) to avoid tiny fragments.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 140:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a monoacylglycerol"

    # Helper: Check if an oxygen atom is an ester oxygen.
    # Here we define ester oxygen as one that (a) is bonded to a carbon (the acyl carbon) that itself has a double bond to an oxygen.
    def is_ester_oxygen(ox, parent_c, mol):
        # For the oxygen (ox) attached to parent carbon, look at its neighbors (excluding the parent)
        for nbr in ox.GetNeighbors():
            if nbr.GetIdx() == parent_c.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                # Look for a double bond from this neighbor (the acyl carbon) to some oxygen.
                for ac_nbr in nbr.GetNeighbors():
                    if ac_nbr.GetAtomicNum() == 8 and ac_nbr.GetIdx() != ox.GetIdx():
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), ac_nbr.GetIdx())
                        if bond and bond.GetBondTypeAsDouble() == 2.0:
                            return True
        return False

    # Helper: Check if an oxygen atom is a free hydroxyl.
    # We require the oxygen to have at least one hydrogen attached.
    def is_free_oh(ox):
        for nbr in ox.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                return True
        return False

    # Now search for candidate glycerol backbones:
    # We look for acyclic paths of length 2 (three connected carbons, in order: A - B - C)
    valid_candidate_found = False
    reason_details = ""
    num_atoms = mol.GetNumAtoms()
    # Iterate over all atoms; candidate atoms must be carbons.
    for atomA in mol.GetAtoms():
        if atomA.GetAtomicNum() != 6 or atomA.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        idxA = atomA.GetIdx()
        # Look at neighbors of A that are carbons (possible middle atoms)
        for atomB in atomA.GetNeighbors():
            if atomB.GetAtomicNum() != 6 or atomB.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            idxB = atomB.GetIdx()
            # First, require bond A-B NOT be in a ring
            bondAB = mol.GetBondBetweenAtoms(idxA, idxB)
            if bondAB is None or bondAB.IsInRing():
                continue
            # Now, from B, go to a third carbon
            for atomC in atomB.GetNeighbors():
                if atomC.GetIdx() in (idxA, idxB):
                    continue
                if atomC.GetAtomicNum() != 6 or atomC.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                idxC = atomC.GetIdx()
                bondBC = mol.GetBondBetweenAtoms(idxB, idxC)
                if bondBC is None or bondBC.IsInRing():
                    continue
                # We now have a candidate 3-carbon chain: idxA, idxB, idxC.
                # For simplicity, sort indices so that the chain is treated uniquely.
                candidate_chain = [atomA, atomB, atomC]
                
                # For each carbon in the candidate chain, collect oxygen neighbors that are NOT in the chain.
                # We also check that these oxygens are not linked to another backbone carbon.
                properties = {}  # key: atom idx, value: 'ester' or 'free' or None
                for carbon in candidate_chain:
                    cidx = carbon.GetIdx()
                    est_flag = False
                    free_flag = False
                    for nbr in carbon.GetNeighbors():
                        # Skip if neighbor is in candidate chain.
                        if nbr.GetIdx() in [a.GetIdx() for a in candidate_chain]:
                            continue
                        if nbr.GetAtomicNum() == 8:
                            # Determine if this oxygen functions as an ester substituent
                            if is_ester_oxygen(nbr, carbon, mol):
                                est_flag = True
                            elif is_free_oh(nbr):
                                free_flag = True
                    # For a glycerol backbone, a given carbon should be either esterified (acyl substitution) or have a free –OH.
                    # (It is allowed that a carbon could have additional substituents, but we require at least one indicator.)
                    if est_flag and free_flag:
                        # Ambiguous substitution; likely not a clean glycerol backbone.
                        properties[cidx] = "ambiguous"
                    elif est_flag:
                        properties[cidx] = "ester"
                    elif free_flag:
                        properties[cidx] = "free"
                    else:
                        properties[cidx] = None

                # To be a monoacylglycerol, exactly one of the three carbons must be esterified and the other two must be free.
                counts = {"ester": 0, "free": 0, "ambiguous": 0, None: 0}
                for val in properties.values():
                    counts[val] += 1
                if counts["ester"] == 1 and counts["free"] == 2:
                    valid_candidate_found = True
                    reason_details = ("Contains a contiguous acyclic 3‐carbon glycerol backbone in which one carbon bears an "
                                      "acyl (ester) substitution while the other two bear free hydroxyl groups.")
                    break
            if valid_candidate_found:
                break
        if valid_candidate_found:
            break

    if valid_candidate_found:
        return True, reason_details
    else:
        # If no candidate chain was found, give details.
        return False, ("No contiguous acyclic 3‐carbon chain with the expected pattern (exactly one ester substitution and "
                       "two free hydroxyl groups) was found. It may be that the glycerol backbone is absent or modified.")

# For testing (this section will run if the file is executed directly)
if __name__ == "__main__":
    # Some test SMILES (from the provided examples)
    test_smiles_list = [
        "O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # MG(0:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0)
        "CCCCCCCC(=O)OC[C@@H](O)CO",  # 1-octanoyl-sn-glycerol
        "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)C(CO)CO",  # MG(0:0/24:6(...)/0:0)
        "O=C1OC(/C(=C/C(=O)OC[C@H](O)CO)/C)=CC(=C1)OC",  # Pestalotiopyrone I (false positive previously)
        "O(C(=O)CCCCCCCCCCCCCCCCCCCCC)C(CO)CO",  # MG(0:0/22:0/0:0) example
    ]
    for s in test_smiles_list:
        result, reason = is_monoacylglycerol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")