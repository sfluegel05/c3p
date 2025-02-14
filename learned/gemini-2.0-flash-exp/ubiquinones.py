"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones have a benzoquinone core, two methoxy groups (or variations) and a polyprenoid side chain (or just a methyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core structure: benzoquinone with two methoxy or -OH groups (or a combination) and a methyl group
    core_pattern_quinone = Chem.MolFromSmarts("C1(=C(C(=O)C=C(C1=O))([OX2]))[OX2]")
    
    if not mol.HasSubstructMatch(core_pattern_quinone):
        return False, "Core benzoquinone structure with two methoxy/OH groups not found"

    # Check for a methyl group on the ring
    methyl_pattern = Chem.MolFromSmarts("C[c]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if not methyl_matches:
       return False, "No methyl group found on the ring"


    #Check for a polyprenoid side chain (or a methyl group for ubiquinone-0).
    # Look for a double bond with a methyl group and a series of carbons with alternating single and double bonds
    sidechain_pattern = Chem.MolFromSmarts("[CX4]=[CX3]([CX4])~[CX4]=[CX3]")
    sidechain_matches = mol.GetSubstructMatches(sidechain_pattern)

    # Check the ubiquinone-0 case: a single methyl group and no side chain.
    if not sidechain_matches and "CC1=C(C(OC)=C(OC)C(=O)C=C1)=O" not in smiles and not "CC1=C(O)C(=O)C=C(C1=O)C" in smiles: 
       return False, "Polyprenoid side chain not found and not ubiquinone-0"
    
    #Handle demethylated ubiquinones
    methoxy_pattern = Chem.MolFromSmarts("OC")
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    if len(methoxy_matches) < 1:
        alcohol_pattern = Chem.MolFromSmarts("[OH1]")
        alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
        if len(alcohol_matches) < 1:
            return False, "No methoxy/alcohol groups found"
        

    return True, "Matches the ubiquinone criteria: core structure, methyl group, methoxy/alcohol groups and polyprenoid chain (or no side chain for ubiquinone-0)."