"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    
    A tannin is a polyphenolic compound, often with astringent properties, recognized by the presence of complex hydroxylated aromatic rings, various linked phenolic structures, glycoside linkages, and ether/ester functionalities.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Known polyphenolic and glycosidic patterns
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    pyrogallol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")
    galloyl_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(OC(=O)c2cc(O)c(O)c(O)c2)c1")
    ellagic_pattern = Chem.MolFromSmarts("O=C1Oc2c(O)c(O)cc(O)c2C(=O)c2c1c(O)c(O)cc2")
    
    # Additional patterns
    condensed_tannin_pattern = Chem.MolFromSmarts("c1[cH]c(O)c(O)c(O)c1-c2c(O)cc(O)c(O)c2")
    glycoside_pattern = Chem.MolFromSmarts("O[C@@H]1[C@@H](O)[C@@H](O)[C@@H](O)[C@H](O)C1")
    
    # Check for important polyphenolic and glycosidic structures
    phenolic_matches = (
        mol.HasSubstructMatch(catechol_pattern) or
        mol.HasSubstructMatch(pyrogallol_pattern) or
        mol.HasSubstructMatch(galloyl_pattern) or
        mol.HasSubstructMatch(ellagic_pattern) or
        mol.HasSubstructMatch(condensed_tannin_pattern)
    )
    
    if not phenolic_matches:
        return False, "No recognized polyphenolic structures found"

    # Evaluate aromatic ring complexity
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 2:
        return False, "Insufficient aromatic ring complexity for tannin identification"
    
    # Assess molecular size for reasonable complexity
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight possibly too low, but could still be tannin-related"
    
    # Check for ester/ether linkages or glycosidic linkages
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    
    ester_or_ether_matches = mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(ether_pattern) or mol.HasSubstructMatch(glycoside_pattern)
    
    if not ester_or_ether_matches:
        return False, "No ester, ether, or glycosidic linkages found"

    return True, "Contains complex polyphenolic structures with sufficient aromaticity and plausible glycosidic/ester features indicative of tannins"