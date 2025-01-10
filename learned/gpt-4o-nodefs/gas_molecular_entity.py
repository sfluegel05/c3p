"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    This is a heuristic approach considering size and structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_weight = Descriptors.ExactMolWt(mol)
    
    # Heuristic rules:
    # 1. Simple atoms and small molecules typically gaseous
    # 2. Noble gases and small hydrocarbons, halides, perfluorinations, etc.
    
    # Check for noble gases or single atom gases
    single_atom_gases = {"[He]", "[Ne]", "[Ar]", "[Kr]", "[Xe]", "[Rn]"}
    if smiles in single_atom_gases:
        return True, "Single atom noble gas"

    # Check if molecular weight is small (heuristic for gases)
    if mol_weight < 100:
        return True, f"Molecular weight is {mol_weight}, typical for gases"

    # Specific small molecules known to be gases
    small_gas_smarts = [
        Chem.MolFromSmarts("O=C=O"), # CO2
        Chem.MolFromSmarts("N"),     # Ammonia
        Chem.MolFromSmarts("[H]Cl"), # Hydrogen chloride
        Chem.MolFromSmarts("[H]I"),  # Hydrogen iodide
        Chem.MolFromSmarts("CC"),    # Ethane
        Chem.MolFromSmarts("C=C"),   # Ethene
        Chem.MolFromSmarts("C#N")    # Hydrogen cyanide
    ]

    for gas_pattern in small_gas_smarts:
        if mol.HasSubstructMatch(gas_pattern):
            return True, "Matches known small gas molecule pattern"

    # Fallback
    return False, "Does not match gas molecular entity heuristics"