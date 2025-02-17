"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: Endocannabinoids – A class of cannabinoids present in mammalian biological fluids
and tissues that activate cannabinoid receptors.
This classifier uses improved heuristics:
  • Excludes molecules with phosphorus.
  • Requires that the molecule have exactly one acyl group:
       – For N‐acylethanolamines: exactly one amide (C(=O)N) and exactly one ethanolamide fragment (C(=O)NCCO)
         with exactly one nitrogen.
       – For monoacylglycerols/glyceryl ethers: exactly one ester (C(=O)O) and exactly one glycerol fragment (C(CO)CO)
         with no nitrogen.
  • In addition, the overall molecule must have at least 18 carbons, at least 15 rotatable bonds,
    and a molecular weight between 250 and 900 Da.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines whether a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids should contain either:
      • an N‐acylethanolamine head group (with one amide bond, pattern C(=O)NCCO, and exactly one nitrogen), or
      • a monoacylglycerol/glyceryl ether head group (with one ester bond, pattern C(=O)O, and no nitrogen).
    In addition, phosphorus atoms are forbidden; the molecule must be “lipid‐like”
    with at least 18 carbons, at least 15 rotatable bonds, and a molecular weight between 250 and 900 Da.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an endocannabinoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Reject molecules containing phosphorus (atomic num 15)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Contains phosphorus; likely a phospholipid rather than an endocannabinoid"
    
    # Define SMARTS patterns for groups we want to check.
    # For ethanolamides (N‐acylethanolamines) our head group is defined by C(=O)NCCO.
    ethanolamide_pattern = Chem.MolFromSmarts("C(=O)NCCO")
    # A more general amide pattern (to count acyl connections in ethanolamides)
    amide_group = Chem.MolFromSmarts("C(=O)N")
    # For glycerol derivatives the head group fragment is defined by C(CO)CO.
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    # General ester pattern (to count the acyl group in monoacylglycerols)
    ester_group = Chem.MolFromSmarts("C(=O)O")
    
    # Count overall nitrogen atoms.
    nN = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    head_valid = False
    head_reason = ""
    
    # Check for ethanolamide head group.
    eth_matches = mol.GetSubstructMatches(ethanolamide_pattern)
    if eth_matches:
        # Expect exactly one ethanolamide fragment and exactly one nitrogen in the whole molecule.
        if len(eth_matches) == 1 and nN == 1:
            # Count amide groups: should be exactly one.
            n_amide = len(mol.GetSubstructMatches(amide_group))
            if n_amide == 1:
                head_valid = True
                head_reason = "Ethanolamide head group found with correct amide and nitrogen counts (1 each)"
            else:
                head_reason = f"Ethanolamide fragment found but number of amide groups is {n_amide} (expected 1)"
        else:
            head_reason = f"Ethanolamide fragment count = {len(eth_matches)} and nitrogen count = {nN} (expected 1 and 1, respectively)"
    
    # If no valid ethanolamide detected, try glycerol head.
    if not head_valid:
        gly_matches = mol.GetSubstructMatches(glycerol_pattern)
        if gly_matches:
            # For a glycerol-based endocannabinoid, expect exactly one glycerol fragment and no nitrogen.
            if len(gly_matches) == 1 and nN == 0:
                n_ester = len(mol.GetSubstructMatches(ester_group))
                if n_ester == 1:
                    head_valid = True
                    head_reason = "Glycerol head group found with correct ester and nitrogen counts (1 and 0, respectively)"
                else:
                    head_reason = f"Glycerol fragment found but number of ester groups is {n_ester} (expected 1)"
            else:
                head_reason = f"Glycerol fragment count = {len(gly_matches)} and nitrogen count = {nN} (expected 1 and 0, respectively)"
    
    if not head_valid:
        return False, "No valid head group found. " + head_reason

    # 2. Check that the molecule overall has a long fatty acyl chain.
    # Require at least 18 carbons.
    nC = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if nC < 18:
        return False, f"Too few carbon atoms ({nC}); expected at least 18 for a fatty acyl chain."
    
    # 3. Check for chain flexibility: require at least 15 rotatable bonds.
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds < 15:
        return False, f"Too few rotatable bonds ({rot_bonds}); expected at least 15 for a flexible lipid chain."
    
    # 4. Ensure molecular weight is within the typical endocannabinoid range (250 to 900 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for an endocannabinoid."
    if mol_wt > 900:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too high for an endocannabinoid."
    
    # Build the detailed reason.
    reason = (f"Valid head group ({head_reason}), no phosphorus, exactly one acyl linkage, "
              f"{nC} carbons, {rot_bonds} rotatable bonds and a molecular weight of {mol_wt:.1f} Da "
              "indicate a lipid‐like structure typical of an endocannabinoid.")
    
    return True, reason