"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: cannabinoids
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids typically contain notable features such as a resorcinol core,
    long carbon chains with unsaturations, ester/amide bonds, and may include cyclic structures.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cannabinoid features:
    
    # Resorcinol core typically seen in THC, CBD derivatives
    resorcinol_core = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    if mol.HasSubstructMatch(resorcinol_core):
        return True, "Classified as cannabinoid by presence of resorcinol core"
    
    # Amide bound to ethanolamine - similar to endogenous cannabinoids like anandamide
    ethanolamine_amide = Chem.MolFromSmarts("C(=O)NCCO") 
    if mol.HasSubstructMatch(ethanolamine_amide): 
        return True, "Classified as cannabinoid by presence of an amide bound to ethanolamine"
    
    # Long unsaturated hydrocarbon chain connected by ester or ether (common in lipophilic cannabinoids)
    long_chain_with_ester_ether = Chem.MolFromSmarts("C=C(CC)C=C(CC)C=C(CC)C=C(OCC)C=O")
    if mol.HasSubstructMatch(long_chain_with_ester_ether):
        return True, "Classified as cannabinoid by presence of long unsaturated chain with ester or ether"
    
    # Bicyclic group adjacent and bridged to hydrocarbon chain (e.g., some synthetic cannabinoids)
    bicyclic_carbocycle_bridge = Chem.MolFromSmarts("c1ccc2c(c1)Cc3ccccc3C2")
    if mol.HasSubstructMatch(bicyclic_carbocycle_bridge):
        return True, "Classified as cannabinoid by presence of bicyclic derivatives"

    return False, "Lacks distinct structural features typical of cannabinoids"