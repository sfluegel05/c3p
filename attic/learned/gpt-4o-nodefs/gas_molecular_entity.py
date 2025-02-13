"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    This is a heuristic approach considering both size, structure, and certain known gaseous patterns.

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
    
    # Noble gases including isotopes
    noble_gases = {"[He]", "[Ne]", "[Ar]", "[Kr]", "[Xe]", "[Rn]",
                   "[220Rn]", "[219Rn]", "[222Rn]"}
    if smiles in noble_gases:
        return (True, "Single atom noble gas or specific isotope")

    # Check if molecular weight fits within typical gas range
    if mol_weight < 100:
        # Known gaseous patterns (e.g., small hydrocarbons, diatomic gases)
        specific_gas_smarts = [
            Chem.MolFromSmarts("O=C=O"),  # CO2
            Chem.MolFromSmarts("[H][H]"),  # H2
            Chem.MolFromSmarts("[O][O]"),  # O2
            Chem.MolFromSmarts("N"),  # N2
            Chem.MolFromSmarts("CC"),  # Ethane
            Chem.MolFromSmarts("C=C"),  # Ethylene
            Chem.MolFromSmarts("[H]C#[O+]")  # Hydrogen cyanide
        ]
        
        for gas_pattern in specific_gas_smarts:
            if mol.HasSubstructMatch(gas_pattern):
                return (True, "Matches specific small gas molecular pattern")

        return (True, f"Molecular weight is {mol_weight}, typical for gases")
    
    # Further checks could include complexity or functionality, omitted for brevity
    
    return (False, "Does not match gas molecular entity criteria")