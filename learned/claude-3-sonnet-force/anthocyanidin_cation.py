"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: CHEBI:35827 anthocyanidin cation
Definition: Any organic cation that is an aglycon of anthocyanin cation; they are oxygenated derivatives of flavylium (2-phenylchromenylium).
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for positive charge
    if sum(atom.GetFormalCharge() for atom in mol.GetAtoms()) != 1:
        return False, "Molecule does not have a positive charge"
    
    # Check for flavylium core
    flavylium_pattern = Chem.MolFromSmarts("[o+]1c2ccccc2oc1")
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "No flavylium core found"
    
    # Check for common anthocyanidin substitution patterns
    pelargonidin_pattern = Chem.MolFromSmarts("c1ccccc1")
    cyanidin_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")
    delphinidin_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    petunidin_pattern = Chem.MolFromSmarts("c1cc(O)c(OC)cc1")
    malvidin_pattern = Chem.MolFromSmarts("c1cc(OC)c(OC)cc1")
    if not any([mol.HasSubstructMatch(pattern) for pattern in [pelargonidin_pattern, cyanidin_pattern, delphinidin_pattern, petunidin_pattern, malvidin_pattern]]):
        return False, "No common anthocyanidin substitution pattern found"
    
    # Check for presence of glycosidic groups
    glycosidic_pattern = Chem.MolFromSmarts("OC")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    
    # Check for additional substituents or counter ions
    substituent_pattern = Chem.MolFromSmarts("[!#1;!#6;!#8;!#7]")
    substituent_matches = mol.GetSubstructMatches(substituent_pattern)
    
    # Classify based on findings
    if glycosidic_matches and substituent_matches:
        return True, "Contains flavylium core with anthocyanidin substitution pattern, glycosidic groups, and additional substituents"
    elif glycosidic_matches:
        return True, "Contains flavylium core with anthocyanidin substitution pattern and glycosidic groups"
    elif substituent_matches:
        return True, "Contains flavylium core with anthocyanidin substitution pattern and additional substituents"
    else:
        return True, "Contains flavylium core with anthocyanidin substitution pattern"