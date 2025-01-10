"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer with two or more units of hydroxyflavans linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Redefine a more general hydroxyflavan structural pattern
    hydroxyflavan_pattern = Chem.MolFromSmarts("c1(cc(O)ccc1)[C@H]2COc3cc(O)ccc3[C@H]2O")
    if not mol.HasSubstructMatch(hydroxyflavan_pattern):
        return False, "No hydroxyflavan units found"

    # Check for typical proanthocyanidin linkage patterns: 4->8, 4->6
    # Considering common flavonoid oligomer linkages
    linkage_4_8_pattern = Chem.MolFromSmarts("c1(cccc(c1))C(CO)Oc2ccccc2")
    linkage_4_6_pattern = Chem.MolFromSmarts("c1ccccc1[C@H](C(CO)O)c2ccccc2")
    
    if not (mol.HasSubstructMatch(linkage_4_8_pattern) or mol.HasSubstructMatch(linkage_4_6_pattern)):
        return False, "Missing flavan linkages (4->8 or 4->6)"

    # Check for oligomeric structure
    # Counting presence of multiple units as a simple method for oligomer check
    flavan_unit_matches = mol.GetSubstructMatches(hydroxyflavan_pattern)
    if len(favan_unit_matches) < 2:
        return False, f"Only {len(favan_unit_matches)} flavan unit(s) found, need at least 2"

    # (Optional) Check existence of gallate, but not necessary for basic proanthocyanidin
    # gallate_pattern = Chem.MolFromSmarts("OC(=O)c1cc(O)c(O)c(O)c1")  # Example of common modification
    # if any(mol.HasSubstructMatch(patt) for patt in [gallate_pattern]):
    #     print("Contains gallate ester groups")

    return True, "Contains multiple hydroxyflavan units with linkages typical of proanthocyanidins"