"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: CHEBI:51673 prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for benzoquinone or naphthoquinone core
    benzoquinone_pattern = Chem.MolFromSmarts("[#6]1=,:[#6](~[#6]=,:[#6](~[#6]1~[#8])~[#8])~[#8]~[#6]")
    naphthoquinone_pattern = Chem.MolFromSmarts("[#6]1=,:[#6](~[#6]=,:[#6](~[#6]2~[#6]=,:[#6]=,:[#6]=,:[#6]=,:[#6]=,:[#6]2~[#8])~[#8])~[#8]")
    if not mol.HasSubstructMatch(benzoquinone_pattern) and not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No benzoquinone or naphthoquinone core found"
    
    # Look for prenyl side chains
    prenyl_pattern = Chem.MolFromSmarts("[#6]-[#6]=,:[#6]-[#6](-[#6])(-[#6])-[#6](-[#6]=,:[#6]-[#6](-[#6])(-[#6])-[#6])*")
    if not mol.HasSubstructMatch(prenyl_pattern):
        return False, "No prenyl side chains found"
    
    # Check for common functional groups and ring systems
    hydroxy_pattern = Chem.MolFromSmarts("[OX1H]")
    methoxy_pattern = Chem.MolFromSmarts("[OX2C]")
    cyclohexene_pattern = Chem.MolFromSmarts("[#6]1=[#6][#6]=[#6][#6]=[#6][#6]1")
    has_hydroxy = mol.HasSubstructMatch(hydroxy_pattern)
    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)
    has_cyclohexene = mol.HasSubstructMatch(cyclohexene_pattern)
    
    if has_hydroxy or has_methoxy or has_cyclohexene:
        return True, "Contains a quinone core, prenyl side chains, and common functional groups or ring systems for prenylquinones"
    else:
        return True, "Contains a quinone core and prenyl side chains"