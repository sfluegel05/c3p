"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: Endocannabinoids – Cannabinoids present in mammalian tissues and fluids that activate cannabinoid receptors.
Typical endocannabinoids include N‐acylethanolamines (e.g. anandamide, palmitoylethanolamide) and monoacylglycerols (e.g. 2‐arachidonoylglycerol, 2‐arachidonyl glyceryl ether).
This improved classifier uses additional heuristics:
  • It excludes molecules containing phosphorus atoms.
  • It rejects molecules with three or more ester groups.
  • It requires a specific head group pattern and also checks that the overall nitrogen content is exactly what one expects for that motif:
       – for N‐acylethanolamines it requires the ethanolamide pattern (C(=O)NCCO) with exactly one nitrogen atom,
       – for monoacylglycerols/glyceryl ethers it requires the glycerol pattern (C(CO)CO) with no nitrogen atoms.
  • It insists that the whole molecule has a robust fatty acyl chain: a minimum carbon count, enough rotatable bonds and a minimum molecular weight.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines whether a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids typically have either:
      - an N-acylethanolamine head group (e.g. –C(=O)NCCO), with exactly one nitrogen, or 
      - a glycerol head group (e.g. fragment C(CO)CO found in monoacylglycerols/glyceryl ethers), with no nitrogen.
    In addition, the molecule must not contain phosphorus and should have only at most two ester groups.
    Finally, we require that the molecule overall is “lipid-like” (enough carbons, a long flexible chain, and a minimum MW).
    
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
    
    # 1. Exclude molecules containing phosphorus (avoid phospholipids)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Contains phosphorus; likely a phospholipid rather than an endocannabinoid"
    
    # 2. Exclude molecules with 3 or more ester groups (e.g. triglycerides)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) >= 3:
        return False, f"Found {len(ester_matches)} ester groups; likely a triglyceride rather than a monoacyl derivative"
    
    # 3. Check for the typical head group.
    # Ethanolamide head group for N‐acylethanolamines:
    ethanolamide_pattern = Chem.MolFromSmarts("C(=O)NCCO")
    # Glycerol head group for monoacylglycerols or glyceryl ethers:
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    has_ethanolamide = mol.HasSubstructMatch(ethanolamide_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    
    # Count overall heteroatoms
    nN = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # Set a flag if the head group is valid.
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
    
    # 4. Check that the fatty acyl chain is long enough.
    # Count the total number of carbon atoms.
    nC = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if nC < 15:
        return False, f"Too few carbons ({nC}); a long fatty acyl chain is expected"
    
    # 5. Check for chain flexibility: require a minimum number of rotatable bonds.
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds < 4:
        return False, f"Too few rotatable bonds ({rot_bonds}); the fatty chain may be too short or rigid for an endocannabinoid"
    
    # 6. Minimum molecular weight is required (~250 Da lower cutoff).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is lower than expected for an endocannabinoid"
    
    reason = (f"Valid head group ({head_reason}), no phosphorus, fewer than 3 ester groups, " +
              f"{nC} carbons, {rot_bonds} rotatable bonds and a molecular weight of {mol_wt:.1f} Da "
              "indicate a lipid-like structure typical of an endocannabinoid.")
    return True, reason