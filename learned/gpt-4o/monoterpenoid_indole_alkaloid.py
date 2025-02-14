"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid (MIA) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is an MIA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for indole substructure: it includes a benzene fused to a pyrrole
    indole_pattern = Chem.MolFromSmarts('c1c[cH]cc2c1[nH]cc2')
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "Indole moiety not found"

    # Check for additional MIA complex scaffold characteristics - partial feature checks
    # (Sekologanin-like terpenoid complexity)
    is_large = Descriptors.MolWt(mol) > 300  # MIAs are typically larger
    has_oxygen = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())  # oxygens present
    has_multiple_chirals = sum(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for atom in mol.GetAtoms()) > 2
    
    # Deduce if structure complexity hints at a terpene-like organization
    if not (is_large and has_oxygen and has_multiple_chirals):
        return False, "Does not meet the expected structure complexity of MIAs"
    
    return True, "Contains indole moiety and exhibits complexity indicative of monoterpenoid indole alkaloids"

# Example usage:
smiles = "C/C=C\1/CN2CC[C@]34C5=CC(C=6C=C7C(=CC6OC)[C@]89CCN%10C/C(=C\C)/[C@](CCC8%10N7C)([C@]9(C(=O)OC)[H])[H])=C(C=C5N(C)[C@]4(C2C[C@@]1([C@@]3(C(=O)OC)[H])[H])[H])OC"
print(is_monoterpenoid_indole_alkaloid(smiles))