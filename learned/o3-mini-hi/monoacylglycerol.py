"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: Monoacylglycerol
Definition: A glyceride in which any one of the three hydroxyl groups (on a contiguous, acyclic three–carbon glycerol backbone)
is esterified (i.e. bears an –O–C(=O)–R group) while the other two bear free –OH.
Heuristic approach (improved):
  1. Parse the molecule and add explicit hydrogens.
  2. Require a minimal molecular weight (here >140 Da).
  3. Find candidate glycerol backbones by looking for a central sp³ carbon (not in a ring)
     that is connected to two distinct sp³ carbon neighbors (also acyclic) that are not bonded to each other.
  4. For the candidate chain (terminal1, central, terminal2), for each carbon check the oxygen substituents
     (neighbors that are not part of the candidate chain). Classify an oxygen as “ester” if it is bonded to a carbon that has a double bond to oxygen.
     Classify an oxygen as “free OH” if it carries at least one explicit hydrogen.
  5. Accept a candidate if exactly one of the three carbons exhibits an ester substitution and the other two each have free hydroxyl groups.
  
Note: This method is heuristic and may mis‐classify molecules with more complex substitution patterns.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    A monoacylglycerol must contain a contiguous acyclic 3–carbon glycerol backbone
    (primary–secondary–primary) in which exactly one position is esterified (i.e. bears an acyl substitution)
    and the other two positions are free (i.e. bear an –OH group).
    
    Args:
        smiles (str): SMILES string.
        
    Returns:
        bool: True if the molecule is classified as monoacylglycerol, False otherwise.
        str: Reason for the decision.
    """
    # Parse SMILES and add explicit hydrogens so that –OH groups are clear.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Guard: require a minimal molecular weight (e.g. >140 Da) to avoid fragments.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 140:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a monoacylglycerol"
    
    # Helper: Return True if this oxygen (ox) acts as part of an ester substituent on the attached carbon.
    def is_ester_oxygen(ox, attached_carbon):
        # Check oxygen (ox) neighbors (except for the attached carbon). 
        # An ester oxygen should be attached to a carbon (the acyl carbon) that in turn is double-bonded to an oxygen.
        for nbr in ox.GetNeighbors():
            if nbr.GetIdx() == attached_carbon.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:  # carbon
                bond = mol.GetBondBetweenAtoms(ox.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                # Look at all neighbors of this potential acyl carbon looking for a double-bonded oxygen.
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetIdx() == ox.GetIdx():
                        continue
                    if subnbr.GetAtomicNum() == 8:
                        bn = mol.GetBondBetweenAtoms(nbr.GetIdx(), subnbr.GetIdx())
                        if bn is not None and bn.GetBondType() == Chem.BondType.DOUBLE:
                            return True
        return False

    # Helper: Return True if oxygen carries at least one hydrogen (free hydroxyl).
    def is_free_oh(ox):
        for nbr in ox.GetNeighbors():
            if nbr.GetAtomicNum() == 1:  # hydrogen
                return True
        return False

    # We now search for candidate glycerol backbones.
    # We require the backbone to be acyclic. In a “true” glycerol the central carbon is connected to two terminal carbons,
    # and these two terminals must not be bonded to each other.
    candidate_found = False
    reason_details = ""
    # Iterate over candidate central carbons.
    for central in mol.GetAtoms():
        if central.GetAtomicNum() != 6:
            continue
        if central.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        # Skip if the central atom is in a ring.
        if central.IsInRing():
            continue
        # Get neighboring carbons that are sp3 and not in ring
        term_neighbors = []
        for nbr in central.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and not nbr.IsInRing():
                term_neighbors.append(nbr)
        # We require at least two distinct terminal carbons.
        if len(term_neighbors) < 2:
            continue
        # Now consider unique pairs among the terminal neighbors.
        n = len(term_neighbors)
        for i in range(n):
            for j in range(i+1, n):
                term1 = term_neighbors[i]
                term2 = term_neighbors[j]
                # Ensure terminal carbons are not directly bonded (this enforces the branching pattern)
                if mol.GetBondBetweenAtoms(term1.GetIdx(), term2.GetIdx()) is not None:
                    continue
                # Our candidate glycerol backbone: [term1, central, term2]
                candidate = [term1, central, term2]
                # For each backbone carbon, check substituents that are oxygens (but not those that are part of the backbone).
                classification = {}  # key: carbon idx, value: "ester", "free", "ambiguous" or None
                for carbon in candidate:
                    poss_types = {"ester": False, "free": False}
                    for nbr in carbon.GetNeighbors():
                        # Exclude backbone atoms.
                        if nbr.GetIdx() in [a.GetIdx() for a in candidate]:
                            continue
                        if nbr.GetAtomicNum() == 8:  # oxygen
                            if is_ester_oxygen(nbr, carbon):
                                poss_types["ester"] = True
                            elif is_free_oh(nbr):
                                poss_types["free"] = True
                    if poss_types["ester"] and poss_types["free"]:
                        classification[carbon.GetIdx()] = "ambiguous"
                    elif poss_types["ester"]:
                        classification[carbon.GetIdx()] = "ester"
                    elif poss_types["free"]:
                        classification[carbon.GetIdx()] = "free"
                    else:
                        classification[carbon.GetIdx()] = None
                # For a valid monoacylglycerol, exactly one must be ester and the other two free.
                num_ester = list(classification.values()).count("ester")
                num_free  = list(classification.values()).count("free")
                if num_ester == 1 and num_free == 2:
                    candidate_found = True
                    reason_details = ("Contains a contiguous acyclic 3‐carbon glycerol backbone (primary–secondary–primary) "
                                      "in which one carbon bears an acyl (ester) substitution and the other two bear free –OH groups.")
                    break
            if candidate_found:
                break
        if candidate_found:
            break

    if candidate_found:
        return True, reason_details
    else:
        return False, ("No contiguous acyclic 3‐carbon chain with the expected pattern (exactly one ester substitution and "
                       "two free hydroxyl groups) was found. It may be that the glycerol backbone is absent or modified.")

# For testing purposes (block executed if module is run directly).
if __name__ == "__main__":
    # A list of test SMILES (some examples from the provided outcomes)
    test_smiles_list = [
        # True positives:
        "O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # MG(0:0/22:6(...)/0:0)
        "CCCCCCCC(=O)OC[C@@H](O)CO",  # 1-octanoyl-sn-glycerol
        "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)C(CO)CO",  # MG(0:0/24:6(...)/0:0)
        "O=C1OC(/C(=C/C(=O)OC[C@H](O)CO)/C)=CC(=C1)OC",  # Pestalotiopyrone I (should be false)
        # Some false negatives (expected not to classify as monoacylglycerol):
        "O(C(=O)CCCCCCCCCCCCCCCCCCCCC)C(CO)CO",  # MG(0:0/22:0/0:0)
    ]
    for s in test_smiles_list:
        res, reason = is_monoacylglycerol(s)
        print(f"SMILES: {s}\nClassification: {res}\nReason: {reason}\n")