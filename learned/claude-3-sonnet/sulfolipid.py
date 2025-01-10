"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: sulfolipid
A compound containing a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Neutralize charges for consistent processing
    uncharger = Chem.rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)

    # Look for sulfonic acid group (C-SO3H) - must have C-S bond
    sulfonic_pattern = Chem.MolFromSmarts("[CX4][S](=[OX1])(=[OX1])[OX2H,OX1-]")
    sulfate_pattern = Chem.MolFromSmarts("[OX2][S](=[OX1])(=[OX1])[OX2H,OX1-]")
    
    if not (mol.HasSubstructMatch(sulfonic_pattern) or mol.HasSubstructMatch(sulfate_pattern)):
        return False, "No sulfonic/sulfate group found"

    # Count carbons and check for minimum lipid size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:  # Lowered threshold to catch simpler sulfolipids
        return False, "Too few carbons for a lipid structure"

    # Look for long carbon chains (at least 8 carbons)
    long_chain = Chem.MolFromSmarts("[CH2,CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "No long carbon chains found"

    # Calculate molecular weight - lowered threshold
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:  # Lowered threshold to include simpler sulfolipids
        return False, "Molecular weight too low for sulfolipid"

    # Look for characteristic lipid features
    features = []
    
    # Check for ester linkages
    ester_pattern = Chem.MolFromSmarts("[#6]C(=O)O[#6]")
    if mol.HasSubstructMatch(ester_pattern):
        features.append("ester bonds")
    
    # Check for amide linkages
    amide_pattern = Chem.MolFromSmarts("[#6]C(=O)N[#6]")
    if mol.HasSubstructMatch(amide_pattern):
        features.append("amide bonds")
    
    # Check for sphingosine backbone
    sphingosine_pattern = Chem.MolFromSmarts("[CH2]O[CH]([CH]([CH]O)O)[CH]([CH](O)/C=C/[#6])N")
    if mol.HasSubstructMatch(sphingosine_pattern):
        features.append("sphingosine backbone")

    if not features:
        return False, "Missing characteristic lipid features"

    # Build reason string
    reason = "Contains sulfonic/sulfate group with lipid features"
    if features:
        reason += f" including {', '.join(features)}"

    return True, reason