"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:36000 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is an ester where the acyl chain (acid component) is derived from lauric acid (dodecanoic acid, a saturated 12-carbon fatty acid).

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

    # Define lauric acid acyl fragment (12 carbons including carbonyl carbon)
    lauric_acid_smiles = 'CCCCCCCCCCC(=O)'  # 11 carbons + carbonyl
    lauric_acid_mol = Chem.MolFromSmiles(lauric_acid_smiles)

    # Identify ester groups: pattern [C](=O)O[C]
    ester_pattern = Chem.MolFromSmarts('[C:1](=O)[O:2][#6]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    for match in ester_matches:
        carbonyl_c_idx = match[0]
        ester_o_idx = match[1]
        alkyl_c_idx = match[2]

        # Break the bond between the ester oxygen and the alkyl carbon
        bond_to_break = mol.GetBondBetweenAtoms(ester_o_idx, alkyl_c_idx).GetIdx()
        fragmented_mol = Chem.FragmentOnBonds(mol, [bond_to_break], addDummies=False)

        # Get the fragments
        frags = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
        # The acyl fragment should contain the carbonyl carbon
        acyl_frag = None
        for frag in frags:
            atom_indices = [atom.GetIdx() for atom in frag.GetAtoms()]
            if carbonyl_c_idx in atom_indices:
                acyl_frag = frag
                break
        if acyl_frag is None:
            continue  # Could not find acyl fragment, continue to next ester group

        # Compare the acyl fragment to lauric acid acyl group
        if acyl_frag.GetNumAtoms() != lauric_acid_mol.GetNumAtoms():
            continue  # Different number of atoms, can't be lauric acid

        # Use molecular fingerprints for comparison
        acyl_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(acyl_frag, radius=2)
        lauric_acid_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(lauric_acid_mol, radius=2)

        similarity = Chem.DataStructs.FingerprintSimilarity(acyl_fp, lauric_acid_fp)
        if similarity == 1.0:
            return True, "Contains dodecanoate ester group"

    return False, "No dodecanoate ester groups with lauric acid found"