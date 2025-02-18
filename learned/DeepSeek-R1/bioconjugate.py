"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate (contains at least one biological component and a covalent linker).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define valid biological component patterns
    amino_acid = Chem.MolFromSmarts("[NH2,NH3+]-[CH](-[C](=[O])[OH,O-])-*")
    sugar = Chem.MolFromSmarts("[C]1O[C@H](O)[C@H](O)[C@H](O)[C@H]1O")
    nucleobase = Chem.MolFromSmarts("[nH]1cnc2ncnc12")
    coa = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    fatty_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H,O-]")
    
    # Check for biological components (skip None patterns)
    bio_matches = 0
    for pattern in [amino_acid, sugar, nucleobase, coa, fatty_acid]:
        if pattern and mol.HasSubstructMatch(pattern):
            bio_matches += 1
    
    # Define covalent linker patterns
    disulfide = Chem.MolFromSmarts("[S][S]")
    thioether = Chem.MolFromSmarts("[S]([#6])[#6]")
    ester = Chem.MolFromSmarts("[O][CX3](=O)[#6]")
    amide = Chem.MolFromSmarts("[CX3](=O)[NX3H]")
    
    # Check for linkers (skip None patterns)
    linker_found = any(
        pattern and mol.HasSubstructMatch(pattern)
        for pattern in [disulfide, thioether, ester, amide]
    )
    
    # Require at least one biological component and one linker
    if bio_matches >= 1 and linker_found:
        return True, "Contains biological component(s) and covalent linker"
    
    # Special case for CoA-like structures with thioester bonds
    thioester = Chem.MolFromSmarts("[#6][S][CX3](=O)[#6]")
    if thioester and mol.HasSubstructMatch(thioester):
        if coa and mol.HasSubstructMatch(coa):
            return True, "Contains CoA-like structure with thioester bond"
    
    return False, "Does not meet bioconjugate criteria"