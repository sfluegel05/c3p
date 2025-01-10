"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.

    A tannin is a polyphenolic compound, often with astringent properties, recognized by the presence of complex hydroxylated aromatic rings and sometimes glycoside bonds.
    
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
    
    # SMARTS patterns for catechol and pyrogallol (known phenolic structures)
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    pyrogallol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")
    galloyl_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(OC(=O)c2cc(O)c(O)c(O)c2)c1")
    
    # Check for any known phenolic structures
    phenolic_matches = (
        mol.GetSubstructMatches(catechol_pattern) or
        mol.GetSubstructMatches(pyrogallol_pattern) or
        mol.GetSubstructMatches(galloyl_pattern)
    )
    
    if not phenolic_matches:
        return False, "No polyphenolic (catechol or pyrogallol-like) structures found"
    
    # Check molecular complexity - number of rings
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 3:
        return False, "Insufficient aromatic ring complexity for tannin identification"
    
    # Check for large molecular size (proxied by molecular weight)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for typical tannin structure"
    
    # Check for ester/ether linkages which are common in glycosides
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    has_ester_or_ether = mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(ether_pattern)
    
    if not has_ester_or_ether:
        return False, "No ester or ether linkages found, indicating insufficient glycosidic linkage or basal complexity"
    
    return True, "Contains complex polyphenolic structures with sufficient aromaticity and plausible glycosidic/ester features indicative of tannins"

# Additional metadata or characterization could be added if necessary