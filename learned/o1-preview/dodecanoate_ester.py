"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:<id> dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is an ester where the carboxylic acid component is lauric acid (dodecanoic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester functional group SMARTS pattern
    ester_smarts = '[C](=O)O[C]'
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester groups found"

    # Define lauric acid acyl group (dodecanoyl) SMARTS
    lauryl_acyl_smarts = 'CCCCCCCCCCC(=O)'  # 11 carbons + carbonyl
    lauryl_acyl = Chem.MolFromSmarts(lauryl_acyl_smarts)

    # Iterate over ester bonds
    for match in ester_matches:
        carbonyl_c_idx = match[0]
        o_idx = match[1]
        alkoxy_c_idx = match[2]

        # Get the acyl side (carboxylic acid part)
        acyl_fragment = Chem.FragmentOnBonds(mol, [mol.GetBondBetweenAtoms(o_idx, alkoxy_c_idx).GetIdx()], addDummies=True)
        acyl_mol = Chem.GetMolFrags(acyl_fragment, asMols=True)[0]

        # Check if the acyl fragment matches lauric acid acyl group
        if acyl_mol.HasSubstructMatch(lauryl_acyl):
            # Further validate that the acyl chain length is exactly 12 carbons
            c_chain = [atom for atom in acyl_mol.GetAtoms() if atom.GetAtomicNum() == 6]
            if len(c_chain) == 12:
                return True, "Contains lauric acid ester group"
            else:
                continue  # Not the correct chain length

    return False, "No lauric acid ester groups found"