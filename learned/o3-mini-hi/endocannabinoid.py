"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: Endocannabinoids – A class of cannabinoids present in mammalian biological fluids
and tissues that activate cannabinoid receptors.
This classifier uses improved heuristics:
  • Excludes molecules with phosphorus.
  • Requires that the molecule have exactly one acyl linkage:
       – For N‐acylethanolamines: exactly one amide (C(=O)N) with an ethanolamide fragment (C(=O)NCCO) and exactly one nitrogen.
       – For monoacylglycerols/glyceryl ethers: a glycerol fragment that is directly linked (via an ester or an ether) to one long fatty acyl chain.
  • In addition, the overall molecule must be “lipid‐like” with at least 18 carbons, at least 15 rotatable bonds,
    and a molecular weight between 250 and 900 Da.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines whether a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are expected to contain either:
      (a) an N‐acylethanolamine head group defined by an amide bond and an ethanolamine fragment (SMARTS: C(=O)NCCO) 
          with exactly one amide and one nitrogen in the whole molecule, or
      (b) a monoacylglycerol/glyceryl ether head group defined by the presence of a glycerol fragment (SMARTS: OC(CO)CO)
          that is directly linked to a long fatty chain via either an ester bond or an ether linkage.
    In addition, molecules with phosphorus are excluded and the overall molecule’s size must be “lipid‐like”
    (at least 18 carbons, at least 15 rotatable bonds, molecular weight between 250 and 900 Da).

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an endocannabinoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Exclude any molecule containing phosphorus (atomic number 15)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Contains phosphorus; likely a phospholipid rather than an endocannabinoid"

    # --- Define SMARTS patterns ---
    # Ethanolamide head group: should include a C(=O)NCCO fragment.
    ethanolamide_pat = Chem.MolFromSmarts("C(=O)NCCO")
    # A general amide group for counting (C(=O)N)
    amide_pat = Chem.MolFromSmarts("C(=O)N")
    # Glycerol head group (we use a relaxed pattern that ignores chirality): O[C@@H](CO)CO, written without stereo here.
    glycerol_pat = Chem.MolFromSmarts("OC(CO)CO")
    # Ester (acyl) bond pattern: C(=O)O
    ester_pat = Chem.MolFromSmarts("C(=O)O")

    # Count total nitrogen (for ethanolamide, we expect exactly 1)
    nN = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    head_valid = False
    head_type = None  # "ethanolamide" or "glycerol"
    head_reason = ""

    # --- Strategy A: Ethanolamide head group ---
    eth_matches = mol.GetSubstructMatches(ethanolamide_pat)
    if eth_matches:
        # To be an ethanolamide, we expect (ideally) one ethanolamide fragment,
        # exactly one amide group and one nitrogen (the head group nitrogen).
        amide_matches = mol.GetSubstructMatches(amide_pat)
        if len(eth_matches) == 1 and len(amide_matches) == 1 and nN == 1:
            head_valid = True
            head_type = "ethanolamide"
            head_reason = "Ethanolamide head group found with correct amide and nitrogen counts (1 each)"
        else:
            head_reason = (f"Ethanolamide fragment count = {len(eth_matches)}, amide count = {len(amide_matches)}, "
                           f"and nitrogen count = {nN} (expected 1, 1, and 1 respectively)")

    # --- Strategy B: Glycerol-based head group (monoacylglycerol or glyceryl ether) ---
    if not head_valid:
        gly_matches = mol.GetSubstructMatches(glycerol_pat)
        if gly_matches:
            # For glycerol-based endocannabinoids,
            # we require that one of the glycerol fragments is linked to exactly one acyl chain.
            # We'll check two possibilities:
            # Option 1. The acyl chain is attached as an ester:
            ester_matches = mol.GetSubstructMatches(ester_pat)
            acyl_count = 0
            # Count ester bonds that appear to connect a glycerol fragment to an acyl chain.
            for gmatch in gly_matches:
                # gmatch is a tuple of atom indices corresponding to the glycerol fragment.
                gset = set(gmatch)
                for ematch in ester_matches:
                    # In an ester match, the oxygen (the second atom) is the one connecting to the alcohol.
                    # If that oxygen is part of our glycerol fragment, then count that ester as the acyl linkage.
                    if ematch[1] in gset:
                        acyl_count += 1
            if acyl_count == 1 and nN == 0:
                head_valid = True
                head_type = "glycerol (ester)"
                head_reason = "Glycerol head group found with one ester acyl linkage and no nitrogen"
            else:
                # Option 2. If no ester bond was found from the glycerol fragment, try to detect an ether linkage.
                # We check each glycerol fragment: for each oxygen atom in the fragment, consider neighbors
                # not in the fragment. If one such neighbor is carbon and is part of a long alkyl chain, count it.
                def count_ether_acyl(gmatch):
                    count = 0
                    for idx in gmatch:
                        atom = mol.GetAtomWithIdx(idx)
                        if atom.GetSymbol() == "O":
                            for nbr in atom.GetNeighbors():
                                if nbr.GetIdx() not in gmatch and nbr.GetAtomicNum() == 6:
                                    # Count carbons in the branch starting from nbr (simple DFS, not counting rings)
                                    seen = set()
                                    stack = [nbr]
                                    branch_carbons = 0
                                    while stack:
                                        a = stack.pop()
                                        if a.GetIdx() in seen:
                                            continue
                                        seen.add(a.GetIdx())
                                        if a.GetAtomicNum() == 6:
                                            branch_carbons += 1
                                        # Only follow non-oxygen neighbors to stay in the hydrocarbon chain
                                        for nn in a.GetNeighbors():
                                            if nn.GetIdx() not in seen and nn.GetAtomicNum() == 6:
                                                stack.append(nn)
                                    # If the branch has at least 6 carbons, consider it an acyl chain.
                                    if branch_carbons >= 6:
                                        count += 1
                    return count

                ether_acyl_total = 0
                for gmatch in gly_matches:
                    ether_count = count_ether_acyl(gmatch)
                    ether_acyl_total += ether_count
                if ether_acyl_total == 1 and nN == 0:
                    head_valid = True
                    head_type = "glycerol (ether)"
                    head_reason = "Glycerol head group found with one ether-linked acyl chain and no nitrogen"
                else:
                    head_reason = (f"Glycerol head group found but acyl linkage count from ester bonds = {acyl_count} "
                                   f"and ether-linked count = {ether_acyl_total} (expected exactly 1) "
                                   f"with nitrogen count = {nN} (expected 0)")

    if not head_valid:
        return False, "No valid head group found. " + head_reason

    # --- Global lipid criteria ---
    # Total carbon atoms (require at least 18)
    nC = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if nC < 18:
        return False, f"Too few carbon atoms ({nC}); expected at least 18 for a fatty acyl chain."

    # At least 15 rotatable bonds
    n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rot < 15:
        return False, f"Too few rotatable bonds ({n_rot}); expected at least 15 for a flexible lipid chain."

    # Molecular weight in [250, 900] Da
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 250:
        return False, f"Molecular weight ({mw:.1f} Da) is too low for an endocannabinoid."
    if mw > 900:
        return False, f"Molecular weight ({mw:.1f} Da) is too high for an endocannabinoid."

    overall_reason = (f"Valid head group ({head_reason}), no phosphorus, exactly one acyl linkage "
                      f"({head_type} based), {nC} carbons, {n_rot} rotatable bonds and a molecular weight of {mw:.1f} Da "
                      "indicate a lipid‐like structure typical of an endocannabinoid.")
    
    return True, overall_reason