"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide consists of four amino acid residues connected by three peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least 3 peptide bonds (-C(=O)N-)
    peptide_bond_pattern = Chem.MolFromSmarts("[-CX3](=[OX1])N")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bonds) != 3:
         return False, f"Found {len(peptide_bonds)} peptide bonds, need exactly 3"

    # Check for a terminal carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxy_matches) < 1:
        return False, f"Missing terminal carboxylic acid group"

    # Check for a terminal amine group, avoid counting the nitrogen of peptide bond
    amine_group_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]")
    amine_matches = mol.GetSubstructMatches(amine_group_pattern)
    
    #Exclude the nitrogens that form a peptide bond
    valid_amine = 0
    for am in amine_matches:
        is_in_peptide = False
        for pb in peptide_bonds:
            if am[0] in pb:
                is_in_peptide = True
        if not is_in_peptide:
           valid_amine += 1
    if valid_amine < 1:
        return False, f"Missing terminal amine group"
    
    
    # Count the number of amino acids. Each peptide bond will have 2 carbonyl carbons and one nitrogen.
    # There will be n+1 C=O carbons, and n+1 nitrogen atoms, where n is the # of peptide bonds.
    # So we expect 4 C=O and 4 Nitrogen atoms
    carbon_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    carbon_count = len(mol.GetSubstructMatches(carbon_pattern))

    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    if carbon_count != 4:
        return False, f"Incorrect carbon count. Found {carbon_count}, must be 4."
    if nitrogen_count != 4:
        return False, f"Incorrect nitrogen count. Found {nitrogen_count}, must be 4."



    return True, "Contains four amino acid residues connected by three peptide linkages"