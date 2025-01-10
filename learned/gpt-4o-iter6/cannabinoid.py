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
    This function checks for structures typical in cannabinoids such as cyclic
    structures with oxygen, long chains with conjugated double bonds or ether linkages,
    and amides linked to ethanolamines.
    
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
    
    # Pattern for cyclic heteroaromatic rings common in THC, CBD derivatives
    heterocyclic_oxygen = Chem.MolFromSmarts("c1cc(O)c(C)c(O)c1")
    if mol.HasSubstructMatch(heterocyclic_oxygen):
        return True, "Classified as cannabinoid by presence of a heteroaromatic ring with oxygen"
    
    # Amide bound to ethanolamine, indicative of anandamide derivatives
    ethanolamine_amide = Chem.MolFromSmarts("CC(=O)NCCO")
    if mol.HasSubstructMatch(ethanolamine_amide):
        return True, "Classified as cannabinoid by presence of an amide bound to ethanolamine"
    
    # Long unsaturated hydrocarbon chain and ether-like structures
    long_ether_chain = Chem.MolFromSmarts("C=CCCCCCCC(=O)OCC")
    if mol.HasSubstructMatch(long_ether_chain):
        return True, "Classified as cannabinoid by presence of ether-linked unsaturated chains"

    # Include tricyclic and polycyclic patterns found in Delta-9-THC and analogs
    tricyclic_pattern = Chem.MolFromSmarts("C1CCC2=C(C1)C=CC3=C2C=CC=C3")
    if mol.HasSubstructMatch(tricyclic_pattern):
        return True, "Classified as cannabinoid by presence of a polycyclic structure typical of THC"

    # Consider additional cannabinoid-like patterns based on known metabolites
    polyol_chain = Chem.MolFromSmarts("C(O)C(O)CC(=O)OCCO")
    if mol.HasSubstructMatch(polyol_chain):
        return True, "Classified as cannabinoid by presence of polyol chains similar to endocannabinoids"
    
    return False, "Lacks distinct structural features typical of cannabinoids"