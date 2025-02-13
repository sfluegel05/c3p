"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is any organic compound with a boiling point <= 250°C under standard atmospheric pressure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a VOC, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = Descriptors.MolWt(mol)

    # Calculate logP as an indicator of potential volatility
    logp = Crippen.MolLogP(mol)

    # Check the length of the carbon chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check for specific functional groups
    alcohol_pattern = Chem.MolFromSmarts("O")
    ether_pattern = Chem.MolFromSmarts("O-C")
    alkene_pattern = Chem.MolFromSmarts("C=C")
    alkyne_pattern = Chem.MolFromSmarts("C#C")
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I,F]")

    # Adjusted logic based on molecular weight, logP, and functional groups
    if (mol.HasSubstructMatch(alcohol_pattern) or 
        mol.HasSubstructMatch(ether_pattern) or 
        mol.HasSubstructMatch(alkene_pattern) or 
        mol.HasSubstructMatch(alkyne_pattern) or 
        mol.HasSubstructMatch(halogen_pattern)):
        
        # Check molecular weight and logP as heuristics for volatility
        if mol_wt <= 400 and logp < 5:
            return True, "Contains VOC functional group and has acceptable molecular weight and logP"

    # Use carbon chain length criteria with a more lenient range for potential VOCs
    if carbon_count <= 22:  # Increased carbon limit to accommodate larger VOCs
        return True, "Carbon chain length and structure suggest potential volatility"

    return False, "Structure does not indicate VOC features"