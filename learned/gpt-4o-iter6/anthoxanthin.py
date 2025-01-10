"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments with flavone, flavonol, or isoflavone cores
    and various oxygen substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Comprehensive pattern set for anthoxanthin core structures
    flavonoid_core_patterns = [
        Chem.MolFromSmarts("c1cc2oc(=O)cc(c2c1)-c1ccc(-c2ccccc2)o1"),  # Flavone core
        Chem.MolFromSmarts("c1cc2oc(=O)c(c2c1)-c1ccccc1"),             # Flavonol core
        Chem.MolFromSmarts("c1c(O)c2c(c1=O)cc(O)cc2"),                 # Isoflavone core
        Chem.MolFromSmarts("c1cc2c(ccc1)oc(=O)c2-c1ccccc1"),           # General chromone derivatives
    ]
    
    core_match = any(mol.HasSubstructMatch(core) for core in flavonoid_core_patterns)
    if not core_match:
        return False, "Flavonoid core structure not found"

    # Check for diversity of oxygenation
    # Consider combinations of hydroxyl, methoxy, and sugar derivatives
    oxygenated_patterns = [
        Chem.MolFromSmarts("[OH]"),                  # Hydroxyl groups
        Chem.MolFromSmarts("[OX2H]"),                # Methoxy groups
        Chem.MolFromSmarts("[CH2][OH]"),             # Sugars
        Chem.MolFromSmarts("c-O-c"),                 # Ethers in aromatic rings
        Chem.MolFromSmarts("c-O[*]"),                # Oxygens linking sugar moieties
    ]
    oxy_matches = set(mol.GetSubstructMatches(Chem.MolFromSmarts("|".join([Chem.MolToSmarts(pat) for pat in oxygenated_patterns]))))

    if len(oxy_matches) < 3:
        return False, f"Insufficient variety of oxygen substitutions, found {len(oxy_matches)} types. Need at least 3."

    return True, "Contains anthoxanthin characteristics with flavonoid scaffold and diverse oxygen substitutions"