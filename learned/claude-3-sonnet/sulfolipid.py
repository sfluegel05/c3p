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

    # Look for C-S bond connected to sulfonic acid group
    # [#6] represents any carbon, S(=O)(=O)[O,O-] represents sulfonic acid/sulfonate
    sulfonic_cs_pattern = Chem.MolFromSmarts("[#6][S](=[OX1])(=[OX1])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(sulfonic_cs_pattern):
        return False, "No carbon-sulfur bond to sulfonic acid group found"

    # Exclude molecules where the sulfur is only connected to oxygen (sulfates)
    sulfate_pattern = Chem.MolFromSmarts("[OX2][S](=[OX1])(=[OX1])[OX2,OX1-]")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    sulfonic_matches = mol.GetSubstructMatches(sulfonic_cs_pattern)
    if len(sulfate_matches) > len(sulfonic_matches):
        return False, "Contains sulfate ester rather than sulfonic acid C-S bond"

    # Check for lipid characteristics
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, "Too few carbons for a lipid structure"

    # Look for long carbon chains (at least 8 carbons)
    long_chain = Chem.MolFromSmarts("[CH2,CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "No long carbon chains found"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:  # Lowered threshold to catch simpler sulfolipids
        return False, "Molecular weight too low for sulfolipid"

    # Count rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:  # Lowered threshold
        return False, "Too rigid for typical sulfolipid"

    # Look for common structural features (but don't require them)
    sugar_pattern = Chem.MolFromSmarts("[CH1,CH2]1O[CH1][CH1][CH1][CH1]1")
    sphingosine_pattern = Chem.MolFromSmarts("[CH2]O[CH1][CH1][CH1][CH1]C[NH2,NH]")
    ester_pattern = Chem.MolFromSmarts("[#6]C(=O)O[#6]")
    amide_pattern = Chem.MolFromSmarts("[#6]C(=O)N[#6]")
    
    features = []
    if mol.HasSubstructMatch(sugar_pattern):
        features.append("sugar moiety")
    if mol.HasSubstructMatch(sphingosine_pattern):
        features.append("sphingosine backbone")
    if mol.HasSubstructMatch(ester_pattern):
        features.append("ester bonds")
    if mol.HasSubstructMatch(amide_pattern):
        features.append("amide bonds")

    # Build reason string
    reason = "Contains sulfonic acid group with C-S bond and lipid characteristics"
    if features:
        reason += f" including {', '.join(features)}"

    return True, reason