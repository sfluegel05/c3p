"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: Monoacylglycerol
Definition: A glyceride in which any one of the three hydroxyl groups (on a contiguous, acyclic three–carbon glycerol backbone)
             is esterified (i.e. bears an –O–C(=O)–R group) while the other two bear free –OH groups.
Heuristic approach (improved):
  1. Parse the molecule and add explicit hydrogens.
  2. Require a minimal molecular weight (>140 Da) to avoid fragments.
  3. Search for a contiguous, acyclic chain of three sp³ carbons (candidate glycerol backbone)
     where the center carbon is connected to exactly two other sp³ carbons that are not directly bonded.
  4. For each backbone carbon, collect oxygen substituents (neighbors not in the backbone).
     (For “true” glycerol, each backbone carbon should have exactly one such oxygen.)
  5. Classify the oxygen group on each carbon as:
       • ester if the oxygen is connected (besides to the backbone) to a carbon that carries a double bond to another oxygen,
       • free hydroxyl if it carries at least one explicit hydrogen.
  6. Accept a candidate backbone if exactly one of the three positions shows an ester substituent and the other two show free –OH.
Note: This method is heuristic and may mis‐classify molecules with complex substitution.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    A monoacylglycerol must have a contiguous acyclic 3–carbon glycerol backbone
    (a linear chain of three sp3 carbons, with the central carbon bonded to exactly two carbon neighbors that are not bonded to each other)
    in which exactly one position is esterified (i.e. bears an –O–C(=O)–R group) and the other two positions bear free –OH groups.
    
    Args:
        smiles (str): SMILES string.
        
    Returns:
        bool: True if the molecule is classified as a monoacylglycerol, False otherwise.
        str: Reason for the decision.
    """
    # Parse SMILES and add explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Guard: require a minimal molecular weight to avoid fragments
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 140:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a monoacylglycerol"
    
    # Helper: Classify an oxygen substituent attached to a candidate backbone carbon.
    def classify_oxygen(ox, attached_carbon):
        # Check if ox acts as an ester oxygen:
        # It must be connected (aside from the backbone carbon) to a carbon (acyl carbon) that is double-bonded to at least one oxygen.
        for nbr in ox.GetNeighbors():
            if nbr.GetIdx() == attached_carbon.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:  # carbon
                bond = mol.GetBondBetweenAtoms(ox.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                # Look for a double-bonded oxygen on this acyl carbon
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetIdx() == ox.GetIdx():
                        continue
                    if subnbr.GetAtomicNum() == 8:
                        bn = mol.GetBondBetweenAtoms(nbr.GetIdx(), subnbr.GetIdx())
                        if bn is not None and bn.GetBondType() == Chem.BondType.DOUBLE:
                            return "ester"
        # If not ester, check if it is a free hydroxyl (has at least one explicit hydrogen neighbor)
        for nbr in ox.GetNeighbors():
            if nbr.GetAtomicNum() == 1:  # hydrogen
                return "free"
        return "unknown"
    
    # Now search for candidate glycerol backbones:
    # A candidate is a contiguous, acyclic chain of three sp3 carbons (a --C--C--C-- sequence)
    # The center carbon must have at least two carbon neighbors; we require exactly two here, and the two must not be directly bonded.
    candidate_found = False
    reason_details = ""
    # Loop over all atoms as potential middle (central) carbon of a glycerol backbone
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        if atom.IsInRing():
            continue
        # Get all carbon neighbors that are sp3 and not in a ring
        carbon_neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and not nbr.IsInRing():
                carbon_neighbors.append(nbr)
        # For glycerol, the central carbon in the backbone should connect to exactly 2 carbons,
        # but in cases where additional C-neighbors exist we try all unique pairs.
        if len(carbon_neighbors) < 2:
            continue
        n = len(carbon_neighbors)
        for i in range(n):
            for j in range(i+1, n):
                term1 = carbon_neighbors[i]
                term2 = carbon_neighbors[j]
                # Enforce that the two terminal carbons are not directly bonded (so the chain is linear)
                if mol.GetBondBetweenAtoms(term1.GetIdx(), term2.GetIdx()) is not None:
                    continue
                # Candidate backbone: [term1, central, term2]
                backbone = [term1, atom, term2]
                # Ensure the three atoms are sp3 and acyclic (they already satisfy, but extra safety)
                valid_backbone = True
                for b in backbone:
                    if b.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or b.IsInRing():
                        valid_backbone = False
                        break
                if not valid_backbone:
                    continue
                # For each backbone carbon, collect oxygen substituents that are not in the backbone.
                substituent_types = {}  # key: backbone atom idx, value: "ester", "free", or "unknown"
                valid_candidate = True
                for b in backbone:
                    oxygens = [nbr for nbr in b.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [a.GetIdx() for a in backbone]]
                    # In glycerol, each backbone carbon should have exactly one oxygen substituent.
                    if len(oxygens) != 1:
                        valid_candidate = False
                        break
                    sub_type = classify_oxygen(oxygens[0], b)
                    if sub_type == "unknown":
                        valid_candidate = False
                        break
                    substituent_types[b.GetIdx()] = sub_type
                if not valid_candidate:
                    continue
                # Count the ester and free classifications across the three backbone carbons.
                ester_count = sum(1 for t in substituent_types.values() if t == "ester")
                free_count = sum(1 for t in substituent_types.values() if t == "free")
                if ester_count == 1 and free_count == 2:
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
        return False, ("No contiguous acyclic 3‐carbon chain with exactly one ester and two free hydroxyl substitutions was found. "
                       "It may be that the glycerol backbone is absent or modified.")

# For testing purposes (executed if module is run directly).
if __name__ == "__main__":
    test_smiles_list = [
        # True positives (monoacylglycerols)
        "O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # MG(0:0/22:6(...)/0:0)
        "CCCCCCCC(=O)OC[C@@H](O)CO",  # 1-octanoyl-sn-glycerol
        "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)C(CO)CO",  # MG(0:0/24:6(...)/0:0)
        "O(=C)OC[C@@H](O)CO",  # Should be invalid/simple fragment
        # False negatives/others:
        "O(C(=O)CCCCCCCCCCCCCCCCCCCCC)C(CO)CO",  # MG(0:0/22:0/0:0) – likely not monoacylglycerol (no proper glycerol backbone)
    ]
    for s in test_smiles_list:
        res, reason = is_monoacylglycerol(s)
        print(f"SMILES: {s}\nClassification: {res}\nReason: {reason}\n")