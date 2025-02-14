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

    # Core structure: benzoquinone with two methoxy or -OH groups (or a combination) 
    # and any substituion at position 5
    core_pattern_quinone = Chem.MolFromSmarts("[c]1[c]([OX2H])[c]([OX2H])[c](=O)[c][c]1=O")
    
    if not mol.HasSubstructMatch(core_pattern_quinone):
        return False, "Core benzoquinone structure with two methoxy/OH groups not found"

    #Verify that the oxygen substituents are at positions 2 and 3:
    methoxy_pattern = Chem.MolFromSmarts("[c]1[c]([OX2H])[c]([OX2H])[c](=O)[c][c]1=O")
    matches_methoxy = mol.GetSubstructMatches(methoxy_pattern)
    if len(matches_methoxy) < 1:
        return False, "No methoxy/alcohol groups in positions 2 and 3"

    #Check for a polyprenoid side chain (or a methyl group for ubiquinone-0)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # If there are more than 4 rotatable bonds and 10 carbon atoms, it has a significant side chain
    if n_rotatable > 4: 
       c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
       if c_count < 10:
           return False, "Not enough carbons to form a polyprenoid side chain"
    # else it could be ubiquinone-0, so we should not discard the molecule if it has no sidechain
    
    return True, "Matches the ubiquinone criteria: core structure with methoxy/alcohol groups at positions 2 and 3, and polyprenoid chain (or no chain for ubiquinone-0)."