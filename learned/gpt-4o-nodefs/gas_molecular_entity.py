"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    This uses a heuristic approach considering size, structure, and known small gas patterns.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Get molecular weight
    mol_weight = Descriptors.ExactMolWt(mol)
    
    # Noble gases and isotopes
    noble_gases = {
        "[He]", "[Ne]", "[Ar]", "[Kr]", "[Xe]", "[Rn]", 
        "[220Rn]", "[219Rn]", "[222Rn]"
    }
    if smiles in noble_gases:
        return (True, "Single atom noble gas or specific isotope")

    # Small molecular weight gases
    if mol_weight < 75:
        specific_gas_smarts = [
            Chem.MolFromSmarts("O=C=O"),   # CO2
            Chem.MolFromSmarts("[H][H]"),  # H2
            Chem.MolFromSmarts("[O][O]"),  # O2
            Chem.MolFromSmarts("N#N"),     # N2
            Chem.MolFromSmarts("C#C"),     # Acetylene
            Chem.MolFromSmarts("CC"),      # Ethane
            Chem.MolFromSmarts("C=C"),     # Ethene
            Chem.MolFromSmarts("C#N"),     # HCN
            Chem.MolFromSmarts("ClCl"),    # Dichlorine
            Chem.MolFromSmarts("FF"),      # Difluorine
            Chem.MolFromSmarts("BrBr"),    # Dibromine
            Chem.MolFromSmarts("II")       # Diiodine
        ]
        
        for gas_pattern in specific_gas_smarts:
            if mol.HasSubstructMatch(gas_pattern):
                return (True, "Matches known gas molecular pattern")
    
    # Check for small hydrocarbons and halogenated hydrocarbons
    if mol_weight < 100:
        small_hydrocarbon_smarts = [
            Chem.MolFromSmarts("C"),        # Methane, ethane, propane patterns
            Chem.MolFromSmarts("C=C=C"),    # Allenes
            Chem.MolFromSmarts("C=C"),      # Ethenes like styrene derivatives
            Chem.MolFromSmarts("C#C"),      # Acetylenes, propyne
        ]
        
        for small_hc_pattern in small_hydrocarbon_smarts:
            if mol.HasSubstructMatch(small_hc_pattern):
                return (True, "Matches small hydrocarbon gas pattern")

    return (False, "Does not match gas molecular entity criteria")