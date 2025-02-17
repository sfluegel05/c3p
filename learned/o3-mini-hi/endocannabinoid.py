"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: Endocannabinoids – A class of cannabinoids present in mammalian tissues and fluids that activate cannabinoid receptors.
This classifier uses multiple heuristics:
  • Excludes molecules with phosphorus (to avoid phospholipids).
  • Excludes molecules with three or more ester groups (to avoid triglycerides).
  • Checks for a valid head group:
       – For N‐acylethanolamines: must find an ethanolamide pattern (C(=O)NCCO) and the total nitrogen count must equal 1.
       – For monoacylglycerols/glyceryl ethers: must find a glycerol pattern (C(CO)CO) and the overall nitrogen count must be 0.
  • In addition, the overall molecule must have enough “lipid‐like” character,
    by having at least 18 carbons, at least 15 rotatable bonds,
    and a molecular weight between 250 and 900 Da.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines whether a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids typically have either an N-acylethanolamine head group (e.g. –C(=O)NCCO)
    or a glycerol head group (e.g. fragment C(CO)CO as found in monoacylglycerols/glyceryl ethers).
    In addition, the molecule must have no phosphorus, fewer than 3 ester groups,
    at least 18 carbons, at least 15 rotatable bonds and a molecular weight between 250 and 900 Da.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an endocannabinoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Convert SMILES to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Exclude molecules with phosphorus atoms (avoid phospholipids)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Contains phosphorus; likely a phospholipid rather than an endocannabinoid"
    
    # 2. Exclude molecules with 3 or more ester groups (avoid triglycerides)
    # Using a simple ester SMARTS pattern [C](=O)O
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) >= 3:
        return False, f"Found {len(ester_matches)} ester groups; likely not a monoacyl derivative"
    
    # 3. Check for the typical head group.
    # Ethanolamide head group (for N‐acylethanolamines)
    ethanolamide_pattern = Chem.MolFromSmarts("C(=O)NCCO")
    # Glycerol head group (for monoacylglycerols/glyceryl ethers)
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    has_ethanolamide = mol.HasSubstructMatch(ethanolamide_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    
    # Count overall nitrogen atoms
    nN = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    head_valid = False
    head_reason = ""
    
    if has_ethanolamide:
        if nN == 1:
            head_valid = True
            head_reason = "Ethanolamide head group found with correct nitrogen count (1)"
        else:
            head_reason = f"Ethanolamide-like fragment found but nitrogen count is {nN} (expected 1)"
    if (not head_valid) and has_glycerol:
        if nN == 0:
            head_valid = True
            head_reason = "Glycerol head group found with correct nitrogen count (0)"
        else:
            head_reason = f"Glycerol-like fragment found but nitrogen count is {nN} (expected 0)"
    
    if not head_valid:
        return False, "No valid head group found. " + head_reason
    
    # 4. Ensure the molecule has enough fatty acyl chain.
    # Require at least 18 carbon atoms.
    nC = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if nC < 18:
        return False, f"Too few carbon atoms ({nC}); expected at least 18 for a long fatty acyl chain"
    
    # 5. Check for chain flexibility by requiring a minimum number of rotatable bonds.
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds < 15:
        return False, f"Too few rotatable bonds ({rot_bonds}); expected at least 15 for a flexible lipid chain"
    
    # 6. Ensure the molecular weight is in an acceptable range (250 to 900 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is lower than expected for an endocannabinoid"
    if mol_wt > 900:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is higher than expected for an endocannabinoid"
    
    reason = (f"Valid head group ({head_reason}), no phosphorus, fewer than 3 ester groups, "
              f"{nC} carbons, {rot_bonds} rotatable bonds and a molecular weight of {mol_wt:.1f} Da "
              "indicate a lipid-like structure typical of an endocannabinoid.")
    
    return True, reason