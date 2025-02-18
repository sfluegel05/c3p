"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: Endocannabinoids – Cannabinoids present in mammalian tissues and fluids that activate cannabinoid receptors.
Examples include various N-acylethanolamines (e.g., anandamide, palmitoylethanolamide) and monoacylglycerols (e.g., 2-arachidonoylglycerol).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids typically contain either:
      - an N-acylethanolamine head group (e.g. –C(=O)NCCO), or 
      - a glycerol head group (e.g. –C(CO)CO as seen in 2-arachidonoylglycerol)
    attached to a long unsaturated fatty acyl chain, and have a molecular weight
    suggestive of a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an endocannabinoid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define SMARTS patterns:
    # Pattern for N-acylethanolamine head group: carbonyl immediately bound to NCCO
    ethanolamide_pattern = Chem.MolFromSmarts("C(=O)NCCO")
    
    # Pattern for glycerol head group as found in 2-arachidonoylglycerol: C(CO)CO
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    
    # Check for presence of either typical head group:
    has_ethanolamide = mol.HasSubstructMatch(ethanolamide_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    if not (has_ethanolamide or has_glycerol):
        return False, "No typical ethanolamine amide or glycerol head group found"
    
    # Ensure the molecule contains a carbonyl group (C(=O)) which is expected in acyl derivatives.
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found (expected in acyl derivatives)"
    
    # Check for unsaturation (i.e. C=C); many fatty acyl chains are unsaturated.
    unsat_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(unsat_pattern):
        return False, "No C=C unsaturation found (fatty acyl chain expected to be unsaturated)"
    
    # Count carbon atoms – endocannabinoids typically have long aliphatic chains. 
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Too few carbons ({c_count}); expected long fatty acyl chain"
    
    # Check molecular weight – many endocannabinoids are in the lipid range (>300 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is lower than expected for an endocannabinoid"
    
    # If all basic checks are passed, classify as an endocannabinoid.
    reason = "Contains a typical ethanolamine or glycerol head group with an acyl carbonyl, unsaturation, and a long aliphatic chain."
    return True, reason