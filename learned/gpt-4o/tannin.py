"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.

    A tannin is a polyphenolic compound, often with astringent properties, recognized by the presence of complex hydroxylated aromatic rings, various linked phenolic structures, and sometimes glycoside linkages.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for common polyphenolic structures (expanded)
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    pyrogallol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")
    galloyl_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(OC(=O)c2cc(O)c(O)c(O)c2)c1")
    ellagic_pattern = Chem.MolFromSmarts("O=C1Oc2c(O)c(O)cc(O)c2C(=O)c2c1c(O)c(O)cc2")

    # Check for any known phenolic structures
    phenolic_matches = (
        mol.GetSubstructMatches(catechol_pattern) or
        mol.GetSubstructMatches(pyrogallol_pattern) or
        mol.GetSubstructMatches(galloyl_pattern) or
        mol.GetSubstructMatches(ellagic_pattern)
    )
    
    if not phenolic_matches:
        return False, "No recognized polyphenolic structures found"
    
    # Evaluate aromatic ring complexity
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 3:
        return False, "Insufficient aromatic ring complexity for tannin identification"
    
    # Assess molecular size without excluding potential candidates prematurely
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight possibly too low, but could still be tannin-related"
    
    # Check for ester/ether linkages common in tannins
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    ester_or_ether_matches = mol.GetSubstructMatches(ester_pattern) or mol.GetSubstructMatches(ether_pattern)
    
    if not ester_or_ether_matches:
        return False, "No ester or ether linkages found, indicating insufficient glycosidic linkage or basal complexity"
    
    return True, "Contains complex polyphenolic structures with sufficient aromaticity and plausible glycosidic/ester features indicative of tannins"

# Additional metadata or characterization could be added if necessary