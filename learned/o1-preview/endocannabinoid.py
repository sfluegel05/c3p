"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.

    An endocannabinoid typically has:
    - A long-chain polyunsaturated fatty acid (usually arachidonic acid derivative)
    - A head group like ethanolamine or glycerol
    - Linked via an amide (for ethanolamine) or ester bond (for glycerol)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
    except Chem.rdchem.KekulizeException:
        return False, "Error parsing SMILES string"

    # Define the fatty acid chain pattern (approximate arachidonic acid)
    fatty_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCC=CC=CC=CC=CCCCCC")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No arachidonic acid-like fatty acid chain found"

    # Define ethanolamine head group pattern
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    has_ethanolamine = mol.HasSubstructMatch(ethanolamine_pattern)

    # Define glycerol head group pattern
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)

    if not has_ethanolamine and not has_glycerol:
        return False, "No ethanolamine or glycerol head group found"

    # Check for amide linkage (for ethanolamine)
    if has_ethanolamine:
        amide_pattern = Chem.MolFromSmarts("C(=O)NCCO")
        if not mol.HasSubstructMatch(amide_pattern):
            return False, "No amide linkage to ethanolamine head group"
        else:
            return True, "Endocannabinoid with ethanolamine head group and amide linkage found"

    # Check for ester linkage (for glycerol)
    if has_glycerol:
        ester_pattern = Chem.MolFromSmarts("C(=O)OCC(O)CO")
        if not mol.HasSubstructMatch(ester_pattern):
            return False, "No ester linkage to glycerol head group"
        else:
            return True, "Endocannabinoid with glycerol head group and ester linkage found"

    return False, "Molecule does not match endocannabinoid criteria"