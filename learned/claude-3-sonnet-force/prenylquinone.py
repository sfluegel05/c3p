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
    
    # Look for quinone core (benzoquinone, naphthoquinone, or other ring systems)
    quinone_patterns = [
        Chem.MolFromSmarts("[#6]1=,:[#6](~[#6]=,:[#6](~[#6]1~[#8])~[#8])~[#8]~[#6]"),  # benzoquinone
        Chem.MolFromSmarts("[#6]1=,:[#6](~[#6]=,:[#6](~[#6]2~[#6]=,:[#6]=,:[#6]=,:[#6]=,:[#6]=,:[#6]2~[#8])~[#8])~[#8]"),  # naphthoquinone
        Chem.MolFromSmarts("[#6]1=,:[#6](~[#6]=,:[#6](~[#6]2~[#6]=,:[#6]=,:[#6]=,:[#6][#6]=,:[#6]2~[#8])~[#8])~[#8]"),  # menaquinone
        Chem.MolFromSmarts("[#6]1=,:[#6](~[#6]=,:[#6](~[#6]2~[#6]=,:[#6]=,:[#6][#6]=,:[#6][#6]=,:[#6]2~[#8])~[#8])~[#8]")  # phylloquinone
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in quinone_patterns):
        return False, "No quinone core found"
    
    # Look for prenyl side chains
    prenyl_pattern = Chem.MolFromSmarts("[#6]-[#6]=,:[#6]-[#6](-[#6])(-[#6])-[#6](-[#6]=,:[#6]-[#6](-[#6])(-[#6])-[#6])*")
    if not mol.HasSubstructMatch(prenyl_pattern):
        return False, "No prenyl side chains found"
    
    # Check for long aliphatic chains
    aliphatic_chain_pattern = Chem.MolFromSmarts("[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]")
    has_long_chain = bool(mol.GetSubstructMatches(aliphatic_chain_pattern))
    
    # Check molecular weight range (typically > 300 Da for prenylquinones)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for prenylquinone"
    
    if has_long_chain:
        return True, "Contains a quinone core, prenyl side chains, and long aliphatic chains"
    else:
        return True, "Contains a quinone core and prenyl side chains"