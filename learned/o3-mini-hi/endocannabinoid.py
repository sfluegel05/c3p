"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: Endocannabinoids – Cannabinoids present in mammalian tissues and fluids that activate cannabinoid receptors.
Typical endocannabinoids include N‐acylethanolamines (e.g. anandamide, palmitoylethanolamide) and monoacylglycerols (e.g. 2‐arachidonoylglycerol, 2‐arachidonyl glyceryl ether).
This classifier uses more stringent criteria by:
  • Excluding molecules containing phosphorus (to filter many phospholipids)
  • Rejecting molecules with three or more ester groups (to rule out triglycerides)
  • Requiring either a specific ethanolamide head group (C(=O)NCCO) or a glycerol head group (C(CO)CO)
  • Checking that the molecule overall has a long fatty acyl chain (enough carbon atoms and rotatable bonds)
  • Requiring a minimal molecular weight (~250 Da)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines whether a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids typically have either:
      - an N-acylethanolamine head group (e.g. –C(=O)NCCO), or 
      - a glycerol head group (e.g. fragment C(CO)CO found in monoacylglycerols or glyceryl ethers)
    attached to a long fatty acyl chain.
    
    Additional checks include:
      - Exclusion of molecules that contain phosphorus (to avoid phospholipids).
      - Rejection of molecules with multiple ester groups (e.g. triglycerides).
      - Sufficient overall carbon count and rotatable bonds to be consistent with a lipid.
      - A minimum molecular weight to catch even smaller endocannabinoids.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an endocannabinoid, False otherwise.
        str: Reason for the classification decision.
    """
    
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules containing phosphorus (many phospholipids contain P)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Contains phosphorus; likely a phospholipid rather than an endocannabinoid"
    
    # Exclude molecules with ≥3 ester groups (triglycerides often have three acyl chains)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) >= 3:
        return False, f"Found {len(ester_matches)} ester groups; likely a triglyceride rather than a monoacyl derivative"
    
    # Define head group patterns:
    # Ethanolamide pattern for N-acylethanolamines
    ethanolamide_pattern = Chem.MolFromSmarts("C(=O)NCCO")
    # Glycerol head group found in monoacylglycerols/glyceryl ethers
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    
    has_ethanolamide = mol.HasSubstructMatch(ethanolamide_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    if not (has_ethanolamide or has_glycerol):
        return False, "No typical ethanolamide (C(=O)NCCO) or glycerol (C(CO)CO) head group found"
    
    # Check for a long fatty acyl chain:
    # We require at least 15 carbon atoms in the molecule.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Too few carbons ({c_count}); a long fatty acyl chain is expected"
    
    # Check that the molecule has a flexible chain
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rotatable_bonds < 3:
        return False, f"Too few rotatable bonds ({rotatable_bonds}); the acyl chain may be too short or rigid"
    
    # Check molecular weight - many endocannabinoids are in the lipid range; lower threshold set to catch palmitoyl ethanolamide (~297 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is lower than expected for an endocannabinoid"
    
    # (Optional) We could check for unsaturation, but because some molecules like palmitoyl ethanolamide are saturated,
    # we relax this requirement.
    
    reason = ("Contains a typical ethanolamide or glycerol head group, no extra phosphorus or multiple ester groups, " +
              "and has a long fatty acyl chain with adequate molecular weight and flexibility")
    return True, reason