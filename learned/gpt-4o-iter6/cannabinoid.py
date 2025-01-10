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
    Cannabinoids often have features like a resorcinol core, amide groups, long
    unsaturated carbon chains, and diverse cyclic structures. 
    
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
    
    # Resorcinol core typical of THC, CBD derivatives
    resorcinol_core = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    if mol.HasSubstructMatch(resorcinol_core):
        return True, "Classified as cannabinoid by presence of resorcinol core"
    
    # Amide bound to ethanolamine, indicative of anandamide derivatives
    ethanolamine_amide = Chem.MolFromSmarts("C(=O)NCCO")
    if mol.HasSubstructMatch(ethanolamine_amide):
        return True, "Classified as cannabinoid by presence of an amide bound to ethanolamine"
    
    # Long unsaturated hydrocarbon chain, common in many cannabinoids
    long_unsaturated_chain = Chem.MolFromSmarts("C=CCCCCCCC(=O)OCC")
    if mol.HasSubstructMatch(long_unsaturated_chain):
        return True, "Classified as cannabinoid by presence of long unsaturated hydrocarbon chains"

    # Incorporate cyclohexene and pentyl chain patterns, typical in Delta-9-THC
    thc_type_pattern = Chem.MolFromSmarts("C1CC[C@H]2[C@H](C1)C(=C)C2")
    if mol.HasSubstructMatch(thc_type_pattern):
        return True, "Classified as cannabinoid by presence of cyclohexene-related structures"

    # Presence of ether linkage and polycyclic backbones
    ether_polycyclic = Chem.MolFromSmarts("C1COC2=C(C1)C=CC=C2")
    if mol.HasSubstructMatch(ether_polycyclic):
        return True, "Classified as cannabinoid by presence of ether-linked polycyclic structures"
    
    return False, "Lacks distinct structural features typical of cannabinoids"