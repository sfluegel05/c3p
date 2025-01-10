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

    # Look for sulfonic acid group (-SO3H)
    sulfonic_pattern = Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(sulfonic_pattern):
        return False, "No sulfonic acid group found"

    # Check for lipid characteristics - long carbon chains
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:  # Minimum carbon count for a lipid
        return False, "Too few carbons for a lipid structure"

    # Look for long carbon chains
    carbon_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No long carbon chains found"

    # Look for common sugar components (optional but common in examples)
    sugar_pattern = Chem.MolFromSmarts("[CH1,CH2]1O[CH1][CH1][CH1][CH1]1")
    has_sugar = mol.HasSubstructMatch(sugar_pattern)

    # Look for ester or amide bonds (common in lipids)
    ester_pattern = Chem.MolFromSmarts("[#6]C(=O)O[#6]")
    amide_pattern = Chem.MolFromSmarts("[#6]C(=O)N[#6]")
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)

    if not (has_ester or has_amide):
        return False, "No ester or amide bonds found typical of lipids"

    # Calculate molecular weight - sulfolipids are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for sulfolipid"

    # Count rotatable bonds to verify flexibility characteristic of lipids
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Too rigid for typical sulfolipid"

    # Success message varies based on structural features
    reason = "Contains sulfonic acid group and lipid characteristics"
    if has_sugar:
        reason += " with sugar moiety"
    if has_ester and has_amide:
        reason += ", including ester and amide bonds"
    elif has_ester:
        reason += ", including ester bonds"
    elif has_amide:
        reason += ", including amide bonds"

    return True, reason